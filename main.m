addpath('Vessel_Library\Vessel_Library_Utilities')
addpath('Vessel_Classes\')


datadir = '../hbp_image_computing/superann/out_h1';
imgdir = '/users/keerthi/documents/work/htic/hbp/data/special/pngcache37/';
imgno = '757';

imgpaths = dir([imgdir '/*_' imgno '_compressed..png']);
img = imread([imgdir '/' imgpaths(1).name]);
lmap = imread([datadir '/labelmap/' imgno '.png' ]);
mpp = 16; % um per pix

spur_length = 150; % smaller than this are spurs

regionnames = {};
regionnames{71}='SSTE';
regionnames{77}='VZ';
regionnames{78}='SVZ';
regionnames{80}='SP';
regionnames{43}='cerebral peduncle';

pix_1mm = fix(1000/mpp);

for rgnid = [71,77,78,80,43]

    rgn = imclose(lmap == rgnid,strel('disk',11));

    % reimplement Vessel_Library/center_spline_fit.m

    [vessels,dist_max] = thinned_vessel_segments(rgn,150);
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
    enforce_connectivity = true;
    set_edges_2nd_derivative(vessels, rgn, [], smooth_parallel, smooth_perpendicular, enforce_connectivity);
    
    figure,imshow(img)
    title(regionnames{rgnid})
    % regionnames{rgnid}
    hold on
    for vessid = 1:numel(vessels)
        % vessid
        nprofiles = size(vessels(vessid).im_profiles,1);
        plot(vessels(vessid).side1(:,2),vessels(vessid).side1(:,1),'y.-')
        plot(vessels(vessid).side2(:,2),vessels(vessid).side2(:,1),'m.-')
        
        

        for ii = unique(fix(linspace(1,nprofiles,ceil(nprofiles/pix_1mm) )))
            x1 = vessels(vessid).side1(ii,2);
            y1 = vessels(vessid).side1(ii,1);
        
            x2 = vessels(vessid).side2(ii,2);
            y2 = vessels(vessid).side2(ii,1);
            
            plot([x1,x2],[y1,y2],'r-')
            len = round(sqrt((x2-x1)^2+(y2-y1)^2),2);
            text(x2+5,y2,num2str(mpp*len),"color",'k')
        end
    end
    hold off

end
