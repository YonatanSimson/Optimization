%Toy IRLS
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
x_ls = lsqr(AA, B, 1e-12, 5000);
x_ls = reshape(x_ls, [rows, cols, dim]);

%solution using inverse - takes a while
Xest = ((A'*A) + lambda*(L'*L))\(A'*y);
residual_ls = sum(sum((x_ls(:) - Xest(:)).^2));

X = reshape(x_ls, [rows, cols, dim]);
% displayVolumeSliceGUI(X);


%% IRLS on previous solution
lambda = 0.1;
epsilon = 0.0001;%for reweighting matrix
maxIter = 50;
alpha   = 0.01;
tol     = 1e-10;
[x, cost] = IRLS(A, L, y, x_ls(:), epsilon, lambda, maxIter, alpha, tol);
figure(1); plot(1:length(cost), cost); xlabel('iter'); title('Cost')
X = reshape(x, [rows, cols, dim]);
