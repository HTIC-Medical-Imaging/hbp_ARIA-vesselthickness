function profiledatalist = measure_profiledata(lmapdata, msk, profilelines, mpp)

ssd = @(p1,p2) sum((p2-p1).^2,2); % is sq eucl
eucldist = @(p1,p2) sqrt(ssd(p1,p2));

profiledatalist = struct('st',[],'en',[],'len',[],'ang',[]);

if isempty(msk)
    msk = true(size(lmapdata{1}));
end
for k = 1:length(profilelines)
    X = [profilelines(k).p1(2),profilelines(k).p2(2)];
    Y = [profilelines(k).p1(1),profilelines(k).p2(1)];
    normal = profilelines(k).u;
    ang = atan2(normal(2),normal(1)); % using r=x, c=y
    % line(X,Y)
    for m=1:length(lmapdata)
        [cx,cy,lp] = improfile((lmapdata{m}>0) & msk,double(X),double(Y)); 
        if sum(lp)==0
            pdata = struct('st',[],'en',[],'len',[],'ang',[]);
            
        else
            intv = get_intervals(lp);
    %         readings = intv(:,2)-intv(:,1);
            stlist = [cy(intv(:,1)), cx(intv(:,1))];
            endlist = [cy(intv(:,2)), cx(intv(:,2))];
            readings = eucldist(stlist,endlist)*mpp;
            pdata = struct('st',stlist,'en',endlist,'len',readings,'ang',ang);
            
        end

        profiledatalist(k,m)=pdata;
        
    end
    
end
