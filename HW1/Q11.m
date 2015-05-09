%% Q11
clear; close all;
load('Small\A');
load('Small\y');

tol     = 1e-10;
maxIter = 5000;

rows = 19;
cols = 19;
dim  = 19;
[Dx, Dy, Dz] = CreateDerivativeOperators3D(rows, cols, dim);


L = [Dx; Dy; Dz];
B = [y; zeros(size(L, 1), 1)];
rng(123);
n = rows*cols*dim;
x0 = 0.5*rand(n, 1);

%% solution using CG
[xx, yy, zz] = meshgrid(1:19, 1:19, 1:19);
for toThePow = [-3 -5 -7],
    lambda = 10^toThePow;
    AA = [A; sqrt(lambda)*L];
    x_cg = CG_LS1(x0,AA,B,tol,maxIter);
    X = reshape(x_cg, [rows, cols, dim]);
    % displayVolumeSliceGUI(X);

    figure;
    p = patch(isosurface(xx,yy,zz,X,0.5));
    isonormals(xx,yy,zz,X,p)
    set(p,'FaceColor','r')
    set(p,'EdgeColor','r')
    daspect([1,1,1])
    view(3); axis tight
    camlight
    lighting gouraud
    title(['$ \lambda = 10^{' num2str(toThePow) '} $'],'interpreter','latex')
end
