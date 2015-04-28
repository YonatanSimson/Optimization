%
% [Gx] = gradient(X);
% Gx = Dx*X

function [Dx] = CreateDerivativeOperators1D(X_length)

%% formulate the Dx, Dy matrices
%Dy ~ diff(X, 1, 1)
idx = 1;%row index on output matrix
ii = zeros(X_length*2, 1);
jj = zeros(X_length*2, 1);
v  = zeros(X_length*2, 1);



ii(2*idx-1) = idx;
jj(2*idx-1) = idx;
v(2*idx-1)  = -1;

ii(2*idx) = idx;
jj(2*idx) = idx + 1;
v(2*idx)  = 1;


idx = idx + 1;
for m = 2:X_length-1,

    ii(2*idx-1) = idx;
    jj(2*idx-1) = idx - 1;
    v(2*idx-1)  = -0.5;

    ii(2*idx) = idx;
    jj(2*idx) = idx + 1;
    v(2*idx)  = 0.5;


    idx = idx + 1;
end

ii(2*idx-1) = idx;
jj(2*idx-1) = idx - 1;
v(2*idx-1)  = -1;

ii(2*idx) = idx;
jj(2*idx) = idx;
v(2*idx)  = 1;

Dy = sparse(ii, jj, v);


Dx = sparse(ii, jj, v);
