function [lineslist, innerpts, outpts] = match_contours(c_outer,c_inner,DM)

% int => closer to ventricle (interior)
% ext => closer to surface (exterior)

    
    lbl_outer = label_boundindices(c_outer.val, 1);
    
    cls = 2;
    if ~c_outer.concave
        cls=1;
    end
    outpts=c_outer.bpts(lbl_outer==cls,:); 
    
%---

    lbl_inner = label_boundindices(c_inner.val, 1);
    
    cls = 1;
    
    if ~c_inner.concave
        cls=2;
    
    end
    intv = get_intervals(lbl_inner==cls);
    
    innerpts = c_inner.bpts(lbl_inner==cls,:);

    pplist = fit_splines(c_inner.bpts, intv, 91);
    
    lineslist = [];
    for j = 1:length(pplist)
        if ~isempty(pplist(j).pp)
            selpts = pplist(j).smoothedpts(20:10:end,:);
            selnormals = pplist(j).normals(20:10:end,:);
            % quiver(selpts(:,2),selpts(:,1),selnormals(:,2),selnormals(:,1),3);
            ind = sub2ind(size(DM),fix(selpts(:,1)),fix(selpts(:,2)));
            lengths = DM(ind)*1.5;
            % p2 = p1+ul
            otherpts = get_boundedpts([selpts(:,1) + selnormals(:,1).*lengths, ...
                         selpts(:,2) + selnormals(:,2).*lengths],size(DM));
            
            startpts = get_boundedpts([selpts(:,1)-selnormals(:,1)*5,...
                        selpts(:,2)-selnormals(:,2)*5],size(DM));
            for k=1:length(selpts)
                ln_k=struct('p1',startpts(k,:), 'p2',otherpts(k,:),'u',selnormals(k,:));
                if isempty(lineslist)
                    lineslist=ln_k;
                else
                    lineslist(end+1)=ln_k;
                end
            end
                
%             X = [selpts(:,2)';otherpts(:,2)']; % 2n x 1
%             Y = [selpts(:,1)';otherpts(:,1)'];
%             line(X,Y)
        end
    end

%%
function lbl = label_boundindices(val, niter)
    if nargin==1
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
