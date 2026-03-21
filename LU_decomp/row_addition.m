function output = row_addition(A, i, scalar, j)
    out = A;
    out(j, :) = out(j, :) + scalar * out(i, :);
    output = out;
end