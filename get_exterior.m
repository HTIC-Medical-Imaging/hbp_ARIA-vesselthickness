function idx_exterior = get_exterior(L)

bbox = regionprops(L,'BoundingBox');

ar_boxes = [];
for bi=1:length(bbox)
    box_i =bbox(bi).BoundingBox; 
    ar_boxes(bi) = box_i(3)*box_i(4);
end
[~,idx_exterior] = max(ar_boxes);

