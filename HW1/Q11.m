%% Q11
clear; close all;
load('small\A');
load('small\y');

rows = 19;
cols = 19;
dim  = 19;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);

lambda = 1e-5;

L = [Dx; Dy; Dz];
AA = [A; sqrt(lambda)*L];
B = [y; zeros(size(L, 1), 1)];
%Solution using lsqr
x_ls = lsqr(AA, B, 1e-6, 200);
x_ls = reshape(x_ls, [rows, cols, dim]);
%solution using inverse - takes a while
Xest = ((A'*A) + lambda*(L'*L))\(A'*y);

residual_ls = sum(sum((x_ls(:) - Xest(:)).^2));

X = reshape(x_ls, [rows, cols, dim]);
displayVolumeSliceGUI(X);