function mask = get_mask_outside(sz, xywh)

mask = false(sz);

r1 = floor(xywh(2));
c1 = floor(xywh(1));
r2 = ceil(r1+xywh(4));
c2 = ceil(c1+xywh(3));
mask(r1:r2,c1:c2)=1;
