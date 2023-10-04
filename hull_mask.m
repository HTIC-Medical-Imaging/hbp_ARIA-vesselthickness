function mask = hull_mask(B,masksize,isconvex)
if nargin==2
    isconvex = false;
end
if isconvex
    [hullidx,~] = convhull(B);
else
    hullidx = boundary(B(:,1),B(:,2));
end
pts = B(hullidx,:);
mask = poly2mask(pts(:,2),pts(:,1),masksize(1),masksize(2));
