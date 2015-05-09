%Solving small problem. Includes IRLS
clear; close all;
load('Large\A');
load('Large\y');

rows = 49;
cols = 49;
dim  = 49;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);

L = [Dx; Dy; Dz];

%% IRLS
alpha = 1;%as required in the question
epsilon = 1e-5;%for reweighting matrix
maxIter = 5000;
tol     = 1e-10;
[x, cost] = IRLS(A, L, y, epsilon, alpha, maxIter, tol);
figure(1); plot(1:length(cost), cost); xlabel('iter'); title('Cost')
X = reshape(x, [rows, cols, dim]);

%% display results
figure(1)
displayVolumeSliceGUI(X);

[xx, yy, zz] = meshgrid(1:rows, 1:cols, 1:dim);

figure(2);
p = patch(isosurface(xx,yy,zz,X,0.5));
isonormals(xx,yy,zz,X,p)
set(p,'FaceColor','r')
set(p,'EdgeColor','r')
daspect([1,1,1])
view(3); axis tight
camlight
lighting gouraud
