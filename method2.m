
% bs  | brain name
% ''''''''''''''''
% 147 | FB36
% 141 | FB40
% 203 | FB63

biosampleid=203;
imgno = 394;

% structuremasks folder available in gdrive (along with google sheet of section numbers)
% https://drive.google.com/drive/folders/1-rmj_KNIhZlMbnpa3myh9o_1k5nHhByR?usp=sharing

datadir = ['./structuremasks/bs-' num2str(biosampleid)];

imgpaths = dir([datadir '/' num2str(imgno) '_nissl.png']);
img = imread([datadir '/' imgpaths(1).name]);

lmapdir = [datadir '/' num2str(imgno)];
lmaps = dir([lmapdir '/*.png']);

regionnames = jsondecode(fscanf(fopen([lmapdir '/names.json']),'%s'));

mpp = 16; % um per pix
pix_1mm = fix(1000/mpp);

selectedids = {'11580','10508','10515','10522','10529','10536','10542'};

selnames = {};
lmapdata = {};
combined = [];
for ii = 1:length(selectedids)
    id = selectedids{ii};
    selnames{ii} = regionnames.(['x' id]);
    lmapdata{ii} = imread([lmapdir '/' id '.png']);
    if ii==1
        combined = uint8(lmapdata{ii}>0);
    else
        combined = combined + uint8(lmapdata{ii}>0)*ii;
    end
end


rgn_vz = lmapdata{7};
rgn_cp = lmapdata{3};

sz = size(rgn_cp);

disp('initialized')

%% ROI selection

% msgbox('select CP ROI followed by VZ ROI (or escape)','modal')

% dispimg = uint8(rgn_cp>0)*100+uint8(rgn_vz>0)*200;
% figure(1),imshow(dispimg,[])    
% 
% rect_cp = drawrectangle('Label','CP region','LabelAlpha',0.2);
% wait(rect_cp)
% ROI_cp = rect_cp.Position;
% 
% assert(~isempty(ROI_cp))
% 
% rect_vz = drawrectangle('Label','Nested VZ region','labelalpha',0.2,'DrawingArea',ROI_cp);
% wait(rect_vz)
% ROI_vz = rect_vz.Position;
% if isempty(ROI_vz)
%     disp('use cp ROI for vz')
%     ROI_vz = ROI_cp;
% end

%% cut rgn_vz and rgn_cp across at the centroid, along 
% 1. minor axis-pi/2
% 2. minor axis+pi/2


% box_cp = regionprops(single(rgn_cp>0),'BoundingBox').BoundingBox;
% mid_cp = box_cp(1:2)+box_cp(3:4)./2; % c,r
mid_cp = regionprops(single(rgn_cp>0),'Centroid').Centroid;
mid_cp_rc = [mid_cp(2),mid_cp(1)];
ori = regionprops(single(rgn_cp>0),'Orientation').Orientation*pi/180;
msk_sectors = get_sector_mask(sz,mid_cp_rc,[ori-pi/2,ori-pi/4,ori,ori+pi/4]);

%% display curves

rgn_cp_masked = rgn_cp & msk_sectors;
rgn_vz_masked = rgn_vz & msk_sectors;

[B_cp, L_cp,~,~]=bwboundaries(rgn_cp_masked,8);
[B_vz, L_vz,~,~]=bwboundaries(rgn_vz_masked,8);

dispimg = uint8(rgn_cp_masked)*10+uint8(rgn_vz_masked)*50;
fi=figure(1);
% imshow(dispimg,[]) 
imshow(img)

plotcontours(fi,B_cp,'cp')
plotcontours(fi,B_vz,'vz')
legend show

%% pairing by nearest

pairings = find_pairing(B_vz,B_cp);

% display each pair
for ii=1:length(pairings)
    vzpts = B_vz{pairings(ii).vz};
    cppts = B_cp{pairings(ii).cp};

    figure(2);
    imshow(img)
    hold on
    
    plot(vzpts(:,2),vzpts(:,1),'-','linewidth',2)
    plot(cppts(:,2),cppts(:,1),'-','linewidth',2)
    title(num2str(struct2array(pairings(ii))))
    hold off
    pause
end

%% 
% pairing cp-n and vz-n

% encl = enclosing(B_cp{pairings(idx,2)},B_vz{pairings(idx,1)},sz);

% FIXME - try sgn = -1 (default 1)

% selected(2) = struct('cp',[2,-6],'vz',1,'vz_within_cp',false,'sgn',1);
% to show that 3 is a child of 2

%% find side A of VZ and side B of CP (group_points), then
% make normals to vz using spline fitting (match_contours), then
% measure on lmapdata, then
% find the vz and cp point on each profile line, then
% filter the profile lines by angle and span, then
% plot the lines

outputdir = [datadir '/' num2str(imgno) '_marked2'];
mkdir(outputdir)

for idx = 4:4 % 1:length(pairings)

    sel = pairings(idx)
    selmsk_cp = select_in_labelmatrix(L_cp,sel.cp);
    selmsk_vz = select_in_labelmatrix(L_vz,sel.vz);
    
    % slow function ahead
    [val_cp, cen_cp] = group_points(B_cp(abs(sel.cp)),selmsk_cp,'bbox');
    [val_vz, cen_vz] = group_points(B_vz(abs(sel.vz)),selmsk_vz,'bbox');
    
    
    boximg = 0*rgn_cp;
    boximg(:,1)=1;
    boximg(:,end)=1;
    boximg(1,:)=1;
    boximg(end,:)=1;
    DM = bwdist(boximg);
    
    contour_cp = struct('bpts',B_cp(abs(sel.cp)),'val',val_cp');
    
    contour_vz = struct('bpts',B_vz(abs(sel.vz)),'val',val_vz');
    sgn = 1;
    [profilelines,vz_pts, cp_pts] = match_contours(contour_cp,contour_vz, DM, sgn, sz);
    
    % figure(fi),hold on
    figure(2),imshow(dispimg,[])
    hold on
    plot(vz_pts(:,2),vz_pts(:,1),'rx') %,'linewidth',2)
    plot(cen_vz(2),cen_vz(1),'rs')
    
    plot(cp_pts(:,2),cp_pts(:,1),'bx') %,'linewidth',2)
    plot(cen_cp(2),cen_cp(1),'bs')
    hold off
    
    
    profiledatalist = measure_profiledata(lmapdata, [], profilelines,mpp);
    % bs/secno/structure
    % pt1,pt2,len
    [profilelines2, angles,spans] = update_profilelines(profilelines,profiledatalist);
    
    valid = filter_profilelines(angles, spans);
    
    fi = plotprofiles(combined, profilelines2,profiledatalist,valid);
    writeprofilecsv(outputdir,selnames,num2str(idx),profiledatalist)
    % saveas(fi,sprintf('%s/markings-%d.png',outputdir,idx))
    exportgraphics(gca,sprintf('%s/markings-%d.png',outputdir,idx))
    save(sprintf('%s/data-%d.mat',outputdir,idx),"profiledatalist",'profilelines')

end

%%

%%
% % pairing cp-2 and vz-2
% contour_outer2 = struct('bpts',B_cp{2},'val',val_cp{2},'concave',false);
% contour_inner2 = struct('bpts',B_vz{2},'val',val_vz{2},'concave',false);
% 
% [profilelines2, innerpts2, outpts2] = match_contours(contour_outer2,contour_inner2,DM);
% 
% profiledatalist2 = measure_profiledata(lmapdata,profilelines2,mpp);
% 
% fi = plotprofiles(combined, profilelines2, profiledatalist2);
% 
% figure(fi),hold on
% plot(innerpts2(:,2),innerpts2(:,1),'rx') %,'linewidth',2)
% plot(outpts2(:,2),outpts2(:,1),'bx') %,'linewidth',2)
% hold off
% 
% %%
% outputdir = [datadir '/' num2str(imgno) '_marked2'];
% mkdir(outputdir)
% writeprofilecsv(outputdir,selnames,'2',profiledatalist)
% saveas(fi,[outputdir '/markings-2.png'])
% save([outputdir '/data-2.mat'],"profiledatalist2",'profilelines2')

%% -- older attempt - shortest distance from int to ext

% sqr = @(xv) xv.^2;
% 
% li = 1;
% for ii = 1:30:length(extpts)
%     dv = sum([sqr(extpts(ii,1)-intpts(:,1)), sqr(extpts(ii,2)-intpts(:,2))],2);
%     [~,intmin_i] = min(dv);
%     intpt = intpts(intmin_i,:);
%     X = [intpt(2);extpts(ii,2)];
%     Y = [intpt(1);extpts(ii,1)];
%     line(X,Y)
% %     profilelines(li).X=X;
% %     profilelines(li).Y=Y;
%     li=li+1;
% end
% 
% for ii = 1:10:length(intpts)
%     dv = (intpts(ii,1)-extpts(:,1)).^2 + (intpts(ii,2)-extpts(:,2)).^2;
%     [~,extmin_i] = min(dv);
%     extpt = extpts(extmin_i,:);
%     X = [intpts(ii,2);extpt(2)];
%     Y = [intpts(ii,1);extpt(1)];
%     line(X,Y)
% end

% hold off
