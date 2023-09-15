
% bs  | brain name
% ''''''''''''''''
% 147 | FB36
% 141 | FB40
% 203 | FB63

biosampleid=203;
imgno = 190;

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

%%
[B_cp, val_cp, cen_cp] = group_points(rgn_cp,'bbox');
[B_vz, val_vz, cen_vz] = group_points(rgn_vz,'hull');

%% visualize the contours so that they can be paired between cp and vz
% also identify if any vz polygon is not concave to ventricle

fi=figure(1);
imshow(img)
plotcontours(fi,B_cp,'cp');
plotcontours(fi,B_vz,'vz');
legend show

%%
boximg = 0*rgn_cp;
boximg(:,1)=1;
boximg(:,end)=1;
boximg(1,:)=1;
boximg(end,:)=1;
DM = bwdist(boximg);

% pairing cp-1 and vz-1
contour_outer = struct('bpts',B_cp{1},'val',val_cp{1},'concave',true);
contour_inner = struct('bpts',B_vz{1},'val',val_vz{1},'concave',true);

[profilelines,innerpts, outpts] = match_contours(contour_outer,contour_inner,DM);

profiledatalist = measure_profiledata(lmapdata,profilelines,mpp);

fi = plotprofiles(combined, profilelines, profiledatalist);

figure(fi),hold on
plot(innerpts(:,2),innerpts(:,1),'rx') %,'linewidth',2)
plot(outpts(:,2),outpts(:,1),'bx') %,'linewidth',2)
hold off


%%
% pairing cp-2 and vz-2
contour_outer2 = struct('bpts',B_cp{2},'val',val_cp{2},'concave',false);
contour_inner2 = struct('bpts',B_vz{2},'val',val_vz{2},'concave',false);

[profilelines2, innerpts2, outpts2] = match_contours(contour_outer2,contour_inner2,DM);

profiledatalist2 = measure_profiledata(lmapdata,profilelines2,mpp);

fi = plotprofiles(combined, profilelines2, profiledatalist2);

figure(fi),hold on
plot(innerpts2(:,2),innerpts2(:,1),'rx') %,'linewidth',2)
plot(outpts2(:,2),outpts2(:,1),'bx') %,'linewidth',2)
hold off

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

hold off
