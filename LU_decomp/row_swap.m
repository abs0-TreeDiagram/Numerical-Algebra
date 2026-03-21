function output = row_swap(A, i, j)
    out = A;
    out([i, j], :) = out([j, i], :);
    output = out;
end