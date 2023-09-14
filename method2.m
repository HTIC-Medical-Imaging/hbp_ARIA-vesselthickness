
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

for ii = 1:length(selectedids)
    id = selectedids{ii};
    selnames{ii} = regionnames.(['x' id]);
    lmapdata{ii} = imread([lmapdir '/' id '.png']);
end


rgn_vz = lmapdata{7};
rgn_cp = lmapdata{3};

%%
[B_cp, val_cp, cen_cp] = group_points(rgn_cp,'bbox');
[B_vz, val_vz, cen_vz] = group_points(rgn_vz,'hull');

%%
figure(1),imshow(img)
hold on
extpts = [];
for k=1:length(B_cp)
    boundpts = B_cp{k};
    int_k = val_cp{k}(:,2)==0; % label as interior
    int_k2 = imdilate(int_k,ones(31,1)); % assume B_cp is a trace
    ext_k = val_cp{k}(:,1)>3 & ~int_k2; % label as exterior

%     plot(boundpts(int_k2,2),boundpts(int_k2,1),'rx')
    plot(boundpts(ext_k,2),boundpts(ext_k,1),'bx')
    extpts=cat(1,extpts,boundpts(ext_k,:));
    break
end
intpts = [];
for k=1:length(B_vz)
    boundpts = B_vz{k};
    int_k = val_vz{k}(:,2)==0;
    int_k2 = imdilate(int_k,ones(31,1)); % assume B_cp is a trace
    ext_k = val_vz{k}(:,1)>3 & ~int_k2;
    ext_k2 = imdilate(ext_k,ones(31,1));
    int_k3 = int_k2 & ~ext_k2;
    intv = get_intervals(int_k3);
    pplist = fit_splines(boundpts, intv, 51);
    
    for j = 1:length(pplist)
        if ~isempty(pplist(j).pp)
            selpts = pplist(j).smoothedpts(1:10:end,:);
            selnormals = pplist(j).normals(1:10:end,:);
            quiver(selpts(:,2),selpts(:,1),selnormals(:,2),selnormals(:,1),3);
        end
    end
    plot(boundpts(int_k3,2),boundpts(int_k3,1),'rx')
%     plot(boundpts(ext_k2,2),boundpts(ext_k2,1),'bx')
    plot(cen_vz{k}(2),cen_vz{k}(1),'gs')
    intpts=cat(1,intpts,boundpts(int_k,:));
    break
end
hold off

%% 

sqr = @(xv) xv.^2;

li = 1;
for ii = 1:30:length(extpts)
    dv = sum([sqr(extpts(ii,1)-intpts(:,1)), sqr(extpts(ii,2)-intpts(:,2))],2);
    [~,intmin_i] = min(dv);
    intpt = intpts(intmin_i,:);
    X = [intpt(2);extpts(ii,2)];
    Y = [intpt(1);extpts(ii,1)];
    line(X,Y)
%     profilelines(li).X=X;
%     profilelines(li).Y=Y;
    li=li+1;
end

% for ii = 1:10:length(intpts)
%     dv = (intpts(ii,1)-extpts(:,1)).^2 + (intpts(ii,2)-extpts(:,2)).^2;
%     [~,extmin_i] = min(dv);
%     extpt = extpts(extmin_i,:);
%     X = [intpts(ii,2);extpt(2)];
%     Y = [intpts(ii,1);extpt(1)];
%     line(X,Y)
% end

hold off
