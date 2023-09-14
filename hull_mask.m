function mask = hull_mask(B,masksize)
% [convhullidx,~] = convhull(B);
hullidx = boundary(B(:,1),B(:,2));
pts = B(hullidx,:);
mask = poly2mask(pts(:,2),pts(:,1),masksize(1),masksize(2));
