function pairings = find_pairing(B_vz,B_cp)
    % curve level pairings
    % returns: struct('cp',cpid,'vz',vzid)
    
    pairings = struct('cp',[],'vz',[]);
    cplist = [];
    for ii = 1:length(B_vz)
        dist_ii = [];
        vzpts = B_vz{ii};
        if length(vzpts)<20 % too few, skip
            continue
        end
        for jj = 1:length(B_cp)
            dist_ii(jj) = get_dist_allpairs(vzpts,B_cp{jj});
        end
        [~,near] = min(dist_ii);
        % nearcppts = B_cp{near};
        cplist(end+1)=near;
        pairings(ii) = struct('cp',near,'vz',ii); 
    end

    for ii = setdiff(1:length(B_cp),cplist)
        dist_ii = [];
        cppts = B_cp{ii};
        if length(cppts)<20
            continue
        end
        for jj = 1:length(B_vz)
            dist_ii(jj) = get_dist_allpairs(cppts,B_vz{jj});
        end
        [~,near] = min(dist_ii);
        
        pairings(end+1) = struct('cp',ii,'vz',near);
    end

function dv = get_dist_allpairs(pts_vz,pts_cp)
    ssd = @(p1,p2) sum((p2-p1).^2,2); % is sq eucl
    D = [];
    for ii = 1:length(pts_vz)
         D(ii) = min(ssd(pts_vz(ii,:),pts_cp));
    end
    dv = mean(D);
    
