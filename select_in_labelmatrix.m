function out = select_in_labelmatrix(L, c_arr)

if length(c_arr)==1
    out = L==c_arr(1);
else
    out = false(size(L));
    for c = c_arr
        if c < 1
            continue
        end
        % out = bitxor(out,L==c);
        out = out | (L==c); % since L won't have overlapping labels
    end
end

