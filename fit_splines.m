function pplist = fit_splines(boundarypts,intervals, piece_spacing)

% boundarypts: one boundary [R,C], obtained from bwboundaries 
% intervals: [start,end] indices

    if nargin == 2
        piece_spacing = 21;
    end

    % pplist = {};
    for ni = 1:size(intervals,1)
        if intervals(ni,2)-intervals(ni,1)<5
            % pplist{ni}=[];
            continue
        else
            selpts = boundarypts(intervals(ni,1):intervals(ni,2),:);
            [pp,pd, smoothedpts, normals ] = fit_spline(selpts, piece_spacing);
            pplist(ni) = struct('pp',pp,'pd',pd,'smoothedpts',smoothedpts,'normals',normals);
        end
    end


function [pp,pd, smoothedpts, normals] = fit_spline(selpts, piece_spacing)


    pixel_dist = [0; sqrt(sum((diff(selpts, [], 1).^2).'))'];
    x = cumsum(sqrt(pixel_dist));
    n_pieces = round(sum(pixel_dist) / piece_spacing);
    if n_pieces < 1
        n_pieces = 1;
    end

    x_eval = (0:floor(max(x)))';

    b = linspace(x(1), x(end), n_pieces+1);

    pp = spline(b, selpts'/spline(b, eye(length(b)), x(:)'));

    smoothedpts = ppval(pp,x_eval)';

    % Compute spline derivative
    pd = pp;
    pd.order = pp.order - 1;
    pd.coefs = bsxfun(@times, pd.coefs(:, 1:3), [3, 2, 1]);
    % Evaluate the derivative
    der = ppval(pd, x_eval);

    normals = [0 1; -1 0] * der;
    normals = normalize_vectors(normals)';
    
function u = normalize_vectors(v)
% v is a 2 * n array containing n 2-dimensional vectors
% u is a 2 * n array containing n 2-dimensional unit vectors

    u = bsxfun(@rdivide, v, hypot(v(1,:), v(2,:)));
