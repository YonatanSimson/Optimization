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

Dx = sparse(X_cols*X_rows*X_dim, X_cols*X_rows*X_dim);
Dy = sparse(X_cols*X_rows*X_dim, X_cols*X_rows*X_dim);
[Dx2d, Dy2d] = CreateDerivativeOperators(X_rows, X_cols);

size2d = X_cols*X_rows;

for l = 1:X_dim,
    Dx((l-1)*size2d + (1:size2d), (l-1)*size2d + (1:size2d)) = Dx2d; 
    Dy((l-1)*size2d + (1:size2d), (l-1)*size2d + (1:size2d)) = Dy2d; %#ok<SPRIX>
    
end

%Dz ~ diff(X, 1, 3)
idx = 1;
ii = zeros(size2d*X_dim*2, 1);
jj = zeros(size2d*X_dim*2, 1);
v  = zeros(size2d*X_dim*2, 1);

for m = 1:size2d,
    
    ii(2*idx-1) = idx;
    jj(2*idx-1) = idx;
    v(2*idx-1)  = -1;
    
    ii(2*idx) = idx;
    jj(2*idx) = idx + size2d;
    v(2*idx)  = 1;
    idx = idx + 1;
end

for n = 2:X_dim-1,
    for m = 1:size2d,
        ii(2*idx-1) = idx;
        jj(2*idx-1) = idx - size2d;
        v(2*idx-1)  = -0.5;

        ii(2*idx) = idx;
        jj(2*idx) = idx + size2d;
        v(2*idx)  = 0.5;
        idx = idx + 1;
    end
end

for m = 1:size2d,

    ii(2*idx-1) = idx;
    jj(2*idx-1) = idx;
    v(2*idx-1)  = 1;
    
    ii(2*idx) = idx;
    jj(2*idx) = idx - size2d;
    v(2*idx)  = -1;

    idx = idx + 1;
end

Dz = sparse(ii, jj, v);

