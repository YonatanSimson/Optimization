%
% [Gx, Gy] = gradient(X);
% (Dx*X(:))'
% (Dy*X(:))'

function [Dx, Dy] = CreateDerivativeOperators(X_rows, X_cols)

%% formulate the L matrix
%Dy - diff(X, 1, 1)
idx = 1;
Dy = zeros(X_cols*X_rows);
for n = 1:X_cols,
    Dy(idx, idx)     = -1;
    Dy(idx, idx + 1) = 1;
    idx = idx + 1;
    for m = 2:X_rows-1,
        Dy(idx, idx - 1) = -0.5;
        Dy(idx, idx + 1) = 0.5;
        idx = idx + 1;
    end
    Dy(idx, idx - 1) = -1;
    Dy(idx, idx)     = 1;
    idx = idx + 1;
end

%Dx - diff(X, 1, 2)
idx = 1;
Dx = zeros(X_cols*X_rows);
for m = 1:X_rows,
    Dx(idx, idx)          = -1;
    Dx(idx, idx + X_rows) = 1;
    idx = idx + 1;
end

for n = 2:X_cols-1,
    for m = 1:X_rows,
        Dx(idx, idx - X_rows) = -0.5;
        Dx(idx, idx + X_rows) = 0.5;
        idx = idx + 1;
    end
end

for m = 1:X_rows,
    Dx(idx, idx)          = 1;
    Dx(idx, idx - X_rows) = -1;
    idx = idx + 1;
end
