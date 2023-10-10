function [lineslist, vz_pts, cp_pts] = match_contours(contour_cp, contour_vz, DM, sgn,sp_piece_spacing)

% int => closer to ventricle (interior)
% ext => closer to surface (exterior)
    
    if nargin == 4
        sp_piece_spacing = 51;
    end
    shp = size(DM);

    vz_within_cp = enclosing([contour_cp.bpts;contour_cp.cen],contour_vz.bpts, shp) | within_box([contour_cp.bpts;contour_cp.cen],contour_vz.bpts)
    cls_cp = 2;
    cls_vz = 1;
    if ~vz_within_cp
        cls_cp = 1;
        cls_vz = 2;
        sgn = -1; % FIXME: why?
    end

    cp_pts = [];
    for k = 1:length(contour_cp)
        lbl_cp = label_boundindices(contour_cp(k).val, 2);
    
        cp_pts=cat(1,cp_pts,contour_cp(k).bpts(lbl_cp==cls_cp,:)); 
    end
%---

    vz_pts = [];
    lineslist = [];
    for k=1:length(contour_vz)
        lbl_vz = label_boundindices(contour_vz(k).val, 3);
    
        intv = get_intervals(lbl_vz==cls_vz,sp_piece_spacing);
    
        vz_pts = cat(1,vz_pts,contour_vz(k).bpts(lbl_vz==cls_vz,:));

        pplist = fit_splines(contour_vz(k).bpts, intv, sp_piece_spacing);
    
    
        for j = 1:length(pplist)
            if ~isempty(pplist(j).pp)
                selpts = pplist(j).smoothedpts(20:10:end,:);
                selnormals = pplist(j).normals(20:10:end,:);
                % quiver(selpts(:,2),selpts(:,1),selnormals(:,2),selnormals(:,1),3);
                ind = sub2ind(shp,fix(selpts(:,1)),fix(selpts(:,2)));
    
                length_multiplier = 1.5; % FIXME
                lengths = sgn*DM(ind)*length_multiplier;
                
                % p2 = p1+ul
                otherpts = get_boundedpts([selpts(:,1) + selnormals(:,1).*lengths, ...
                             selpts(:,2) + selnormals(:,2).*lengths],shp);
                
                length_multiplier_back = 10;  % FIXME 
                npix_back = sgn*length_multiplier_back;
                startpts = get_boundedpts([selpts(:,1)-selnormals(:,1)*npix_back,...
                            selpts(:,2)-selnormals(:,2)*npix_back],shp);
    
                for k=1:length(selpts)
                    ln_k=struct('p1',startpts(k,:), 'p2',otherpts(k,:),'u',selnormals(k,:));
                    if isempty(lineslist)
                        lineslist=ln_k;
                    else
                        lineslist(end+1)=ln_k; % FIXME: only place with end+1 append
                    end
                end
                    
    %             X = [selpts(:,2)';otherpts(:,2)']; % 2n x 1
    %             Y = [selpts(:,1)';otherpts(:,1)'];
    %             line(X,Y)
            end
        end

    end

%%
function lbl = label_boundindices(val, niter)
    if nargin==1
        niter = 1;
    end
    if niter < 1
        niter = 1;
    end
    lbl_int = val(:,2)==0;
    lbl_ext = val(:,1)>3;
    for ii = 1:niter
        lbl_int = imdilate(lbl_int,ones(31,1));
        lbl_ext = lbl_ext & ~lbl_int;
    end
    lbl = zeros(size(lbl_int),'uint8');
    % 2=> exterior, 1=> interior, 0=>unclassified
    lbl(lbl_int)=1;
    lbl(lbl_ext)=2;

%%
function RC2 = get_boundedpts(RC,siz)
    nr = siz(1);
    nc = siz(2);
    R = RC(:,1);
    C = RC(:,2);
    R2 = min(nr,max(1,R));
    C2 = min(nc,max(1,C));
    RC2 = [R2,C2];

%%
function out = enclosing(pts_cp,pts_vz, shp)
    
    cp_hull = hull_mask(pts_cp,shp,1);
    pts_vz_lin = sub2ind(shp,pts_vz(:,1),pts_vz(:,2));
    v = cp_hull(pts_vz_lin);
    out = sum(v==1) > sum(v==0);

%%
function out = within_box(pts_cp,pts_vz)

    rc_min = min(pts_cp,[],1);
    rmin = rc_min(1);
    cmin = rc_min(2);
    rc_max = max(pts_cp,[],1);
    rmax = rc_max(1);
    cmax = rc_max(2);

    v = (pts_vz(:,1)>rmin) & (pts_vz(:,1) < rmax) & (pts_vz(:,2)>cmin) & (pts_vz(:,2)<cmax);
    out = sum(v==1) > sum(v==0);
