function main(biosampleid, imgno)

addpath('Vessel_Library\Vessel_Library_Utilities')
addpath('Vessel_Classes\')

% biosampleid=147;
% imgno = 409;

pyenv('Version','../envs/hbp_env/Scripts/python.exe');

pyrunfile(sprintf("./geojson_to_mask.py '%d' '%d'",biosampleid, imgno))

datadir = ['./structuremasks/bs-' num2str(biosampleid)];

imgpaths = dir([datadir '/' num2str(imgno) '_nissl.png']);
img = imread([datadir '/' imgpaths(1).name]);

lmapdir = [datadir '/' num2str(imgno)];
lmaps = dir([lmapdir '/*.png']);

% lmap = imread([datadir '/labelmap/' imgno '.png' ]);
mpp = 16; % um per pix

spur_length = 150; % smaller than this are spurs

regionnames = jsondecode(fscanf(fopen([lmapdir '/names.json']),'%s'));
% regionnames = {};
% regionnames{71}='SSTE';
% regionnames{77}='VZ';
% regionnames{78}='SVZ';
% regionnames{80}='SP';
% regionnames{43}='cerebral peduncle';

pix_1mm = fix(1000/mpp);


outputdir = [datadir '/' num2str(imgno) '_marked'];
if  exist(outputdir,'dir')
    rmdir(outputdir,'s')
end
mkdir(outputdir)

selectedids = ['11580','10508','10515','10522','10529','10536','10542'];
% selectedids = [];

for lm = lmaps(3:end)'
    
    ontoid = lm.name(1:end-4);
    if ~isempty(selectedids) && ~contains(selectedids,ontoid)
%     if ~strcmp(ontoid, '10522')
        continue
    end
    rgnname = regionnames.(['x' ontoid]);
    disp([ontoid ' ' rgnname])
    fp = fopen([outputdir '/' rgnname '_lines.txt'],'wt');

    rgn = imclose(imread([lm.folder '/' lm.name])>0,strel('disk',11));

    % reimplement Vessel_Library/center_spline_fit.m
    try
        ii=0;
        while ii<5
            [vessels,dist_max] = thinned_vessel_segments(rgn,150,2,false,true);
            if dist_max > 0
                break
            else
                disp('dist_max=0; retrying')
                ii=ii+1;
            end
        end
    catch er
        disp('skipping')
        continue
    end
    spline_piece_spacing = 10; % spacing between spline pieces
    vessels = spline_centreline(vessels, spline_piece_spacing);
    
    width = ceil(dist_max * 4);
    if mod(width,2)==0
        width = width + 1;
    end

        
    vessels = make_image_profiles(vessels,uint8(rgn)*200,width,'*linear');
        
    % reimplement Vessel_Library/edges_max_gradient.m
    smooth_parallel = 1;
    smooth_perpendicular = 0.1;
    enforce_connectivity = false;
    set_edges_2nd_derivative(vessels, rgn, [], smooth_parallel, smooth_perpendicular, enforce_connectivity);
    
    fi=figure();
    figure(fi),imshow(img)
    title([num2str(imgno) ' ' rgnname '[' ontoid ']'], 'Interpreter','none')
    % regionnames{rgnid}
    hold on
    for vessid = 1:numel(vessels)
        % vessid
        nprofiles = size(vessels(vessid).im_profiles,1);
        plot(vessels(vessid).side1(:,2),vessels(vessid).side1(:,1),'y.-')
        plot(vessels(vessid).side2(:,2),vessels(vessid).side2(:,1),'m.-')
        plot(vessels(vessid).centre(:,2),vessels(vessid).centre(:,1),'g-')
        
        

        for ii = unique(fix(linspace(1,nprofiles,ceil(nprofiles/pix_1mm) )))
            x1 = vessels(vessid).side1(ii,2);
            y1 = vessels(vessid).side1(ii,1);
        
            x2 = vessels(vessid).side2(ii,2);
            y2 = vessels(vessid).side2(ii,1);
            
            plot([x1,x2],[y1,y2],'r-')
            len = round(sqrt((x2-x1)^2+(y2-y1)^2),2);
            text(x2+5,y2,num2str(mpp*len),"color",'k')
            if ~isnan(len)
                fprintf(fp,'%.2f,%.2f;%.2f,%.2f;%.2f\n',x1,y1,x2,y2,mpp*len);
            end
        end
    end
    hold off
    saveas(fi,[outputdir '/' rgnname '.png'])
    fclose(fp);
    close(fi)
end
