function profiledatalist = measure_profiledata(lmapdata,profilelines, mpp)

ssd = @(p1,p2) sum((p2-p1).^2,2);

profiledatalist = [];

for k = 1:length(profilelines)
    X = [profilelines(k).p1(2),profilelines(k).p2(2)];
    Y = [profilelines(k).p1(1),profilelines(k).p2(1)];
    % line(X,Y)
    for m=1:length(lmapdata)
        [cx,cy,lp] = improfile(lmapdata{m}>0,double(X),double(Y)); 
        if sum(lp)==0
            continue
        end
        intv = get_intervals(lp);
%         readings = intv(:,2)-intv(:,1);
        stlist = [cy(intv(:,1)), cx(intv(:,1))];
        endlist = [cy(intv(:,2)), cx(intv(:,2))];
        readings = sqrt(ssd(stlist,endlist))*mpp;
        pdata = struct('st',stlist,'en',endlist,'len',readings);
        if isempty(profiledatalist)
            profiledatalist = pdata;
        else
            profiledatalist(k,m)=pdata;
        end
    end
    
end
