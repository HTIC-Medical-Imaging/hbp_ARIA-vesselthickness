function fi = plotprofiles(img, profilelines, profiledatalist)

fi=figure();
imshow(img,[])
hold on

for k = 1:length(profilelines)
    linestart = profilelines(k).p1;
    
    linestpts = cell2mat({profiledatalist(k,7).st}'); % start of vz
    if isempty(linestpts)
        line([profilelines(k).p1(2),profilelines(k).p2(2)],...
            [profilelines(k).p1(1),profilelines(k).p2(1)],'color','r')
        continue % FIXME - why happening
    end
    v = sum(abs(linestpts-linestart),2);
    [~,lineminidx] = min(v);

    lineminpt = linestpts(lineminidx,:);
    
    lineenpts = cell2mat({profiledatalist(k,3).en}'); % end of cp
    if isempty(lineenpts)
        line([profilelines(k).p1(2),profilelines(k).p2(2)],...
            [profilelines(k).p1(1),profilelines(k).p2(1)],'color','y')
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
            
            len = pdata.len(1);
            
            ll = sum((lineminpt-pdata.en(1,:)).^2,2);
            if mod(k,5)==0 && ll < linelenmax+1
                plot(pdata.en(1,2),pdata.en(1,1),'y+')

                if mod(k,20)==0
                    text(pdata.en(1,2),pdata.en(1,1),sprintf('%.2f',len),'color','w')
                end
            end
            
            
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
    X = [lineminpt(2),linemaxpt(2)];
    Y = [lineminpt(1),linemaxpt(1)];
    line(X,Y)
%     plot(lineminpt(2),lineminpt(1),'ys')
end

hold off
