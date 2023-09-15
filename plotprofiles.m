function fi = plotprofiles(img, profilelines, profiledatalist)

fi=figure(2);
imshow(img,[])
hold on

for k = 1:length(profilelines)
    linestart = profilelines(k).p1;
%     X = [profilelines(k).p1(2),profilelines(k).p2(2)];
%     Y = [profilelines(k).p1(1),profilelines(k).p2(1)];
%     line(X,Y)
    linelenmax = 0;
    linemaxpt = [];
    for m=1:size(profiledatalist,2)
        pdata = profiledatalist(k,m);
        if ~isempty(pdata.len)
            len = pdata.len(1);
            if mod(k,20)==0
                text(pdata.en(1,2),pdata.en(1,1),sprintf('%.2f',len),'color','w')
            end
            if mod(k,5)==0
                plot(pdata.en(1,2),pdata.en(1,1),'y+')
            end
            ll = sum((linestart-pdata.en(1,:)).^2,2);
            if ll>linelenmax
                linemaxpt=pdata.en(1,:);
                linelenmax = ll;
            end
        end
    end
    X = [linestart(2),linemaxpt(2)];
    Y = [linestart(1),linemaxpt(1)];
    line(X,Y)
end

hold off
