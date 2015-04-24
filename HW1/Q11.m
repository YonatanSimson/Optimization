%% Q11
clear; close all;
load('Small\A');
load('Small\y');

rows = 19;
cols = 19;
dim  = 19;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);

lambda = 1e-3;

L = [Dx; Dy; Dz];
AA = [A; sqrt(lambda)*L];
B = [y; zeros(size(L, 1), 1)];
%Solution using lsqr
x_ls = lsqr(AA, B, 1e-10, 5000);
x_ls = reshape(x_ls, [rows, cols, dim]);
%solution using inverse - takes a while
Xest = ((A'*A) + lambda*(L'*L))\(A'*y);

residual_ls = sum(sum((x_ls(:) - Xest(:)).^2));

X = reshape(x_ls, [rows, cols, dim]);
displayVolumeSliceGUI(X);

%% Large
clear; close all;
load('Large\A');
load('Large\y');

rows = 49;
cols = 49;
dim  = 49;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);

lambda = 1e-5;

L = [Dx; Dy; Dz];
AA = [A; sqrt(lambda)*L];
B = [y; zeros(size(L, 1), 1)];
%Solution using lsqr
x_ls = lsqr(AA, B, 1e-10, 5000);
x_ls = reshape(x_ls, [rows, cols, dim]);

X = reshape(x_ls, [rows, cols, dim]);
displayVolumeSliceGUI(X);

X = reshape(x_ls, [rows, cols, dim]);
displayVolumeSliceGUI(X);