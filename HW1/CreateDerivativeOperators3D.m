%
% [Gx, Gy, Gz] = gradient(X);
% Gx = reshape(Dx*X(:), size(X))
% Gy = reshape(Dy*X(:), size(X))
% Gz = reshape(Dz*X(:), size(X))

function [Dx, Dy, Dz] = CreateDerivativeOperators3D(X_rows, X_cols, X_dim)

%% formulate the Dx, Dy matrices
%Dy ~ diff(X, 1, 1)
idx = 1;%row index on output matrix
ii = zeros(X_cols*X_rows*X_dim*2, 1);
jj = zeros(X_cols*X_rows*X_dim*2, 1);
v  = zeros(X_cols*X_rows*X_dim*2, 1);

Dy = sparse(X_cols*X_rows*X_dim, X_cols*X_rows*X_dim);
% Dy = [];
size2d = X_cols*X_rows;
for l = 1:X_dim,
    [Dx2D, Dy2d] = CreateDerivativeOperators(X_rows, X_cols);
    Dy((l-1)*size2d + (1:size2d), (l-1)*size2d + (1:size2d)) = Dy2d;
end
Dy = sparse(ii, jj, v);

%Dx ~ diff(X, 1, 2)
idx = 1;
ii = zeros(X_cols*X_rows*2, 1);
jj = zeros(X_cols*X_rows*2, 1);
v  = zeros(X_cols*X_rows*2, 1);

for m = 1:X_rows,
    
    ii(2*idx-1) = idx;
    jj(2*idx-1) = idx;
    v(2*idx-1)  = -1;
    
    ii(2*idx) = idx;
    jj(2*idx) = idx + X_rows;
    v(2*idx)  = 1;
    idx = idx + 1;
end

for n = 2:X_cols-1,
    for m = 1:X_rows,
        ii(2*idx-1) = idx;
        jj(2*idx-1) = idx - X_rows;
        v(2*idx-1)  = -0.5;

        ii(2*idx) = idx;
        jj(2*idx) = idx + X_rows;
        v(2*idx)  = 0.5;
        idx = idx + 1;
    end
end

for m = 1:X_rows,

    ii(2*idx-1) = idx;
    jj(2*idx-1) = idx;
    v(2*idx-1)  = 1;
    
    ii(2*idx) = idx;
    jj(2*idx) = idx - X_rows;
    v(2*idx)  = -1;

    idx = idx + 1;
end

Dx = sparse(ii, jj, v);

%Dz~ diff(X, 1, 3)
Dz = zeros(size(Dx));