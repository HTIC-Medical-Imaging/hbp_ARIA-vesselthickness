function [profilelines2, angles, spans] = update_profilelines(profilelines,profiledatalist)
% return augmented struct, and array of angles and spans (for filtering)

profilelines2 = struct('p1',[],'p2',[],'u',[],'ang',[],'span',[]);

angles = [];
spans = [];

for k = 1:length(profilelines)

    profilelines2(k).p1 = profilelines(k).p1;
    profilelines2(k).p2 = profilelines(k).p2;
    normal = profilelines(k).u;
    profilelines2(k).u = normal;
    ang = atan2(normal(2),normal(1)); % r=x, c=y
    profilelines2(k).ang = ang;
    profilelines2(k).span = -1;

    angles(k)=ang;
    spans(k)=-1;
    linestpts = cell2mat({profiledatalist(k,7).st}'); % start of vz
    if isempty(linestpts)
        continue
    end

    linestart = profilelines(k).p1;

    v = sum(abs(linestpts-linestart),2);
    [~,lineminidx] = min(v);

    lineminpt = linestpts(lineminidx,:);
    
    lineenpts = cell2mat({profiledatalist(k,3).en}'); % end of cp
    if isempty(lineenpts)
        continue % FIXME : why is this happening
    end
    
    v = sum(abs(lineenpts-linestart),2);
    [~,lineminidx] = min(v);
    
    linemaxpt = lineenpts(lineminidx,:); 

    if isempty(linemaxpt)
        linelenmax = 0;
    else
        linelenmax = sum((lineminpt-linemaxpt).^2,2); % sq eucl
    end
    linemaxpt_srch = [];

    for m=1:size(profiledatalist,2)
        pdata = profiledatalist(k,m);
        
        if ~isempty(pdata.len)
                    
            ll = sum((lineminpt-pdata.en(1,:)).^2,2);
            
            if isempty(linemaxpt)
                
                if ll>linelenmax
                    linemaxpt_srch=pdata.en(1,:);
                    linelenmax = ll;
                end
            end
        end
    end

    if isempty(linemaxpt)
        linemaxpt = linemaxpt_srch;
    end
%     X = [lineminpt(2),linemaxpt(2)];
%     Y = [lineminpt(1),linemaxpt(1)];
    profilelines2(k).p1 = lineminpt;
    profilelines2(k).p2 = linemaxpt;
    span = sqrt(sum((linemaxpt-lineminpt).^2));
    profilelines2(k).span = span;
    spans(k)=span;
end

