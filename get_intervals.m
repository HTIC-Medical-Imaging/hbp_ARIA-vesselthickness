function intv = get_intervals(selarr)
% find intervals in a bit array (0-1 transitions are starts, 1-0 are ends)
% return a list of (start,end)

d1=diff([0;selarr;0],1);
starts = find(d1==1);
ends = find(d1==-1);
intv = [starts,min(ends,length(selarr))];
