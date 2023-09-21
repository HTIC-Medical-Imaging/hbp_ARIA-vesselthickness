
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

%% ROI selection

msgbox('select CP ROI followed by VZ ROI (or escape)','modal')

dispimg = uint8(rgn_cp>0)*100+uint8(rgn_vz>0)*200;
figure(1),imshow(dispimg,[])    

rect_cp = drawrectangle('Label','CP region','LabelAlpha',0.2);
wait(rect_cp)
ROI_cp = rect_cp.Position;

assert(~isempty(ROI_cp))

rect_vz = drawrectangle('Label','Nested VZ region','labelalpha',0.2,'DrawingArea',ROI_cp);
wait(rect_vz)
ROI_vz = rect_vz.Position;
if isempty(ROI_vz)
    disp('use cp ROI for vz')
    ROI_vz = ROI_cp;
end

%%

msk_cp = get_mask_outside(size(rgn_cp),ROI_cp);
rgn_cp_masked = rgn_cp & msk_cp;
msk_vz = get_mask_outside(size(rgn_vz),ROI_vz);
rgn_vz_masked = rgn_vz & msk_vz;


[B_cp, L_cp,~,~]=bwboundaries(rgn_cp_masked,8);
[B_vz, L_vz,~,~]=bwboundaries(rgn_vz_masked,8);

dispimg = uint8(rgn_cp_masked)*10+uint8(rgn_vz_masked)*50;
fi=figure(1);
% imshow(dispimg,[]) 
imshow(img)

plotcontours(fi,B_cp,'cp')
plotcontours(fi,B_vz,'vz')
legend show

%% 
% pairing cp-n and vz-n
selected = struct('cp',[],'vz',[],'vz_within_cp',true,'sgn',1);

selected(1) = struct('cp',1,'vz',1, 'vz_within_cp',true,'sgn',1);


% FIXME - try sgn = -1 (default 1)

selected(2) = struct('cp',[2,-3],'vz',2,'vz_within_cp',false,'sgn',1);
% to show that 3 is a child of 2

%% Match side A of VZ and side B of CP
% depending on Vz_within_cp (usually true -> A is inner, B is outer)

sel = selected(2)
selmsk_cp = select_in_labelmatrix(L_cp,sel.cp);
selmsk_vz = select_in_labelmatrix(L_vz,sel.vz);
[val_cp, cen_cp] = group_points(B_cp(abs(sel.cp)),selmsk_cp,'pts');
[val_vz, cen_vz] = group_points(B_vz(abs(sel.vz)),selmsk_vz,'bbox');

%%
boximg = 0*rgn_cp;
boximg(:,1)=1;
boximg(:,end)=1;
boximg(1,:)=1;
boximg(end,:)=1;
DM = bwdist(boximg);

contour_cp = struct('bpts',B_cp(abs(sel.cp)),'val',val_cp');

contour_vz = struct('bpts',B_vz(abs(sel.vz)),'val',val_vz');

[profilelines,vz_pts, cp_pts] = match_contours(contour_cp,contour_vz, sel.vz_within_cp, DM, sel.sgn);

% figure(fi),hold on
figure(2),imshow(dispimg,[])
hold on
plot(vz_pts(:,2),vz_pts(:,1),'rx') %,'linewidth',2)
plot(cen_vz(:,2),cen_vz(:,1),'rs')

plot(cp_pts(:,2),cp_pts(:,1),'bx') %,'linewidth',2)
hold off

%%
profiledatalist = measure_profiledata(lmapdata, msk_cp, profilelines,mpp);
% bs/secno/structure
% pt1,pt2,len

fi = plotprofiles(combined, profilelines, profiledatalist);


%%
outputdir = [datadir '/' num2str(imgno) '_marked2'];
mkdir(outputdir)
writeprofilecsv(outputdir,selnames,'1',profiledatalist)
saveas(fi,[outputdir '/markings-1.png'])
save([outputdir '/data-1.mat'],"profiledatalist",'profilelines')

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
