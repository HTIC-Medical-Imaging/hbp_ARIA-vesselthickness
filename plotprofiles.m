function fi = plotprofiles(img, profilelines2, profiledatalist,valid)

fi=figure(3);
imshow(img,[])
hold on

for k = 1:length(profilelines2)
    span = profilelines2(k).span;
    if ~valid(k)
        continue
    end
    if span < 0
        continue
    end
    
    for m=1:size(profiledatalist,2)
        pdata = profiledatalist(k,m);
        
        if ~isempty(pdata.len)
            
            len = pdata.len(1);
            
            % ll = sum((lineminpt-pdata.en(1,:)).^2,2);
            if mod(k,5)==0 % && ll < linelenmax+1
                plot(pdata.en(1,2),pdata.en(1,1),'y+')

                if mod(k,10)==0
                    text(pdata.en(1,2),pdata.en(1,1),sprintf('%.2f',len),'color','r')
                end
            end
            
        end
    end

    X = [profilelines2(k).p1(2),profilelines2(k).p2(2)];
    Y = [profilelines2(k).p1(1),profilelines2(k).p2(1)];
    line(X,Y)
   
    plot(profilelines2(k).p1(2),profilelines2(k).p1(1),'ys')
%     text(profilelines2(k).p1(2),profilelines2(k).p1(1),num2str(k),'color','r')
end

hold off
