function [pt_labels,centroids] = group_points(B_k,L_k,centroid_method)
    % return pt_labels (array of [profile lengths, last 5 sum])
    % for interior points, pt_labels{i}(:,2) == 0
    % for exterior points, pt_labels{i}(:,1) > 3

    if nargin == 2
        centroid_method = 'lmap'; % lmap, bbox, hull, pts
    end

    pt_labels = {};
    centroids = {};
    sz = size(L_k);
    
        
    for k = 1:length(B_k) % nested hole boundaries
        
        boundpts = B_k{k};
        if length(boundpts)<100
            % continue % too small region - possibly erroneous
            return
        end
        props = regionprops(L_k,{'Centroid','BoundingBox'});
        if props(1).BoundingBox(3)<10 && props(1).BoundingBox(4)<10
            % continue % too small box
            return
        end
        if strcmp(centroid_method,'lmap')
            centroid_k = props(1).Centroid;
        elseif strcmp(centroid_method,'pts')
            centroid_k = fliplr(mean(boundpts,1));
        elseif strcmp(centroid_method,'bbox')
            bbox = props(1).BoundingBox;
            centroid_k = bbox(1:2)+bbox(3:4)./2;
        elseif strcmp(centroid_method,'hull')
            hullmsk = hull_mask(B_k,size(L_k));
            centroid_k = regionprops(hullmsk,'centroid').Centroid;
        end
        
        cen_x = centroid_k(1);
        cen_y = centroid_k(2);
        centroids{k} = [cen_y,cen_x]; % r,c
        
        
            
%         profileimg = imdilate(bwperim(L==k),strel("square",2));
        profileimg = L_k;
        pt_labels{k}=[];
        
%         visited = false(sz);
        for ii = 1:length(boundpts)

            pt = boundpts(ii,:); % r,c;
%             if visited(pt(1),pt(2))
%                 continue
%             end
            X = [cen_x;pt(2)];
            Y = [cen_y;pt(1)];

            [cx,cy,lp] = improfile(profileimg,X,Y);
%             vidx = sub2ind(sz,fix(cy),fix(cx));

%             visited(vidx)=1;

            nres = sum(lp);
            if length(lp)>5
                nres2 = sum(lp(end-5:end-1));
            else
                nres2 = 0; % FIXME
            end
            % nres2 = sum(lp2);
            pt_labels{k}(ii,:)=[nres, nres2];
        end    

    end
    centroids = cell2mat(centroids');
