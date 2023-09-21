function intv = get_intervals(selarr,min_slicelen)
% find intervals in a bit array (0-1 transitions are starts, 1-0 are ends)
% return a list of (start,end)
if nargin==1
    min_slicelen = 3;
end
d1=diff([0;selarr;0],1);
starts = find(d1==1);
ends = find(d1==-1);
intv = [max(1,starts-1),min(ends-1,length(selarr))];

sel = (intv(:,2)-intv(:,1)) > min_slicelen;

intv = intv(sel,:);