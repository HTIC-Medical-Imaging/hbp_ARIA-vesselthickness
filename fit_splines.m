function pplist = fit_splines(boundarypts,intervals, piece_spacing)

% boundarypts: one boundary [R,C], obtained from bwboundaries 
% intervals: [start,end] indices

    if nargin == 2
        piece_spacing = 21;
    end

    % pplist = {};
    for ni = 1:size(intervals,1)
        if intervals(ni,2)-intervals(ni,1)<piece_spacing
%             pplist(ni).pp=[];
            continue
        else
            selpts = boundarypts(intervals(ni,1):intervals(ni,2),:);
            [pp,pd, smoothedpts, normals ] = fit_spline(selpts, piece_spacing);
            pplist(ni) = struct('pp',pp,'pd',pd,'pts',selpts,'smoothedpts',smoothedpts,'normals',normals);
        end
    end
    if ~exist('pplist','var')
        pplist = [];
    end


function [pp,pd, smoothedpts, normals] = fit_spline(selpts, piece_spacing)


    pixel_dist = [0; sqrt(sum((diff(selpts, [], 1).^2).'))'];
    x = cumsum(sqrt(pixel_dist));
    n_pieces = round(sum(pixel_dist) / piece_spacing);
    if n_pieces < 1
        n_pieces = 1;
    end
    
    x_eval = (0:floor(max(x)))';

    if n_pieces == 1

        deg = min(3, numel(x)-1);
        % Construct Vandermonde matrices
        V = ones(numel(x), deg+1);
        V_eval = ones(numel(x_eval), deg+1);
        for jj = deg:-1:1
            V(:,jj) = x.*V(:,jj+1);
            V_eval(:,jj) = x_eval.*V_eval(:,jj+1);
        end
        % Do fit to get polynomial coefficients - each column gives the
        % coefficients of a polynomial
        pp = V \ selpts;
        % Evaluate polynomial for centre line
        smoothedpts = transpose(pp' * V_eval');
        % Compute polynomial derivate
        pd = bsxfun(@times, pp(1:end-1, :), (deg:-1:1)');
        % Evaluate derivatives
        der = pd' * V_eval(:, 2:end)';
        
        % Note: This code does the same, but more slowly
%         p1 = polyfit(x, y(:,1), 2);
%         p2 = polyfit(x, y(:,2), 2);
%         cent = [polyval(p1, x_eval), polyval(p2, x_eval)];
%         pd1 = polyder(p1);
%         pd2 = polyder(p2);
%         der = [polyval(pd1, x_eval), polyval(pd2, x_eval)];
    else
        
        b = linspace(x(1), x(end), n_pieces+1);
    
        pp = spline(b, selpts'/spline(b, eye(length(b)), x(:)'));

        smoothedpts = ppval(pp,x_eval)';

        % Compute spline derivative
        pd = pp;
        pd.order = pp.order - 1;
        pd.coefs = bsxfun(@times, pd.coefs(:, 1:3), [3, 2, 1]);
        % Evaluate the derivative
        der = ppval(pd, x_eval);

    end

    normals = [0 1; -1 0] * der;
    normals = normalize_vectors(normals)';
    
function u = normalize_vectors(v)
% v is a 2 * n array containing n 2-dimensional vectors
% u is a 2 * n array containing n 2-dimensional unit vectors

    u = bsxfun(@rdivide, v, hypot(v(1,:), v(2,:)));
