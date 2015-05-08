%Solving small problem. Includes IRLS
clear; close all;
load('Small\A');
load('Small\y');

rows = 19;
cols = 19;
dim  = 19;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);

lambda = 1e-3;
L = [Dx; Dy; Dz];

%Solution using lsqr
AA = [A; sqrt(lambda)*L];
B = [y; zeros(size(L, 1), 1)];
x_ls = lsqr(AA, B, 1e-12, 5000);
x_ls = reshape(x_ls, [rows, cols, dim]);

%solution using inverse - takes a while
Xest = ((A'*A) + lambda*(L'*L))\(A'*y);
residual_ls = sum(sum((x_ls(:) - Xest(:)).^2));

%Solution using CG - as learned in class
x0 = 0.5*rand(rows*cols*dim, 1);
tol     = 1e-13;
MaxIt   = 50;
x_cg = CG_LS(x0,AA,B,tol,MaxIt);%Solve_CG(x0,y,A,lambda,L,1e-13,5000);
disp('Amit''s CG solver accuracy for small problem:')
norm(Xest(:)-x_cg)


%% IRLS on previous solution
lambda = 1;%as required in the question
epsilon = 1e-5;%for reweighting matrix
maxIter = 5000;
tol     = 1e-10;
[x, cost] = IRLS(A, L, y, epsilon, lambda, maxIter, tol);
figure(1); plot(1:length(cost), cost); xlabel('iter'); title('Cost')
X = reshape(x, [rows, cols, dim]);

%% display results
figure(1)
displayVolumeSliceGUI(X);

[xx, yy, zz] = meshgrid(1:19, 1:19, 1:19);

figure(2);
p = patch(isosurface(xx,yy,zz,X,0.5));
isonormals(xx,yy,zz,X,p)
set(p,'FaceColor','r')
set(p,'EdgeColor','r')
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud
