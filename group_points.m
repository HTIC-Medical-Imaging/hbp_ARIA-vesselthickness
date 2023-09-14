function [B,pt_labels,centroids] = group_points(mask,centroid_method)
    % return boundaries, and pt_labels (array of [profile lengths, last 5 sum])
    % for interior points, pt_labels{i}(:,2) == 0
    % for exterior points, pt_labels{i}(:,1) > 3
    if nargin == 1
        centroid_method = 'lmap'; % lmap, bbox, hull, pts
    end
    [B, L,n,~]=bwboundaries(mask,8,'noholes');
    pt_labels = {};
    centroids = {};
    for k = 1:n
        boundpts = B{k};
        props = regionprops(L==k,{'Centroid','BoundingBox'});
        if strcmp(centroid_method,'lmap')
            centroid_k = props(1).Centroid;
        elseif strcmp(centroid_method,'pts')
            centroid_k = fliplr(mean(boundpts,1));
        elseif strcmp(centroid_method,'bbox')
            bbox = props(1).BoundingBox;
            centroid_k = bbox(1:2)+bbox(3:4)./2;
        elseif strcmp(centroid_method,'hull')
            hullmsk = hull_mask(B{k},size(L));
            centroid_k = regionprops(hullmsk,'centroid').Centroid;
        end
        centroids{k}=centroid_k;
        cen_x = centroid_k(1);
        cen_y = centroid_k(2);

        
        pt_labels{k}=[];
            
%         profileimg = imdilate(bwperim(L==k),strel("square",2));
        for ii = 1:length(boundpts)

            pt = boundpts(ii,:); % r,c;

            X = [cen_x;pt(2)];
            Y = [cen_y;pt(1)];

            lp = improfile(L==k,X,Y);
%             lp2 = improfile(L==k,X,Y);
            
            nres = sum(lp);
            nres2 = sum(lp(end-5:end-1));
%             nres2 = sum(lp2);
            pt_labels{k}(ii,:)=[nres, nres2];
        end    

    end
