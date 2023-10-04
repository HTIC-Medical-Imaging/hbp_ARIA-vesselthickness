function out = get_sector_mask(sz,cen_rc,orilist)

out = true(sz);

for ori = orilist
    v = [sin(ori),cos(ori)]; % direction of minor axis
    % x = a cos(ori)

    amin = abs((sz(1)-cen_rc(1))/sin(ori));
    for a = 2:amin-2
        pt = fix(cen_rc - a*v);
        if pt(1)==1 || pt(2)==1 || pt(1)>sz(1)-1 || pt(2)>sz(2)-1
            break
        end
        out(pt(1)-1:pt(1)+1,pt(2)-1:pt(2)+1)=0;
        
    end
    
    amax = abs(cen_rc(1)/sin(ori));
    for a = 2:fix(amax)-2
        
        pt = fix(cen_rc + a*v);
        if pt(1)==1 || pt(2)==1 || pt(1)>sz(1)-1 || pt(2)>sz(2)-1
            break
        end
        out(pt(1)-1:pt(1)+1,pt(2)-1:pt(2)+1)=0;
    end
end