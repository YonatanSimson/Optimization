%% Q18 - Surface normals
clear; close all;

load Images.mat
load LightSources.mat
load mozart.mat

I = double(Images);
L = double(LightSources);
%N - number of pictures
[rows, cols, N] = size(I);
Image_n = zeros(rows, cols, 3);

Lpinv = pinv(L);
for n = 1:cols,
    for m = 1:rows,
        Imn = squeeze(I(m, n, :));
        Image_n(m, n, :) = Lpinv*Imn;
    end
end

%find Zx, Zy, N = (-Zx, -Zy, 1)
p = -Image_n(:, :, 1)./(Image_n(:, :, 3));
q = -Image_n(:, :, 2)./(Image_n(:, :, 3));


%% Q19 - Jacobi method
rows_p2 = rows + 2; 
cols_p2 = cols + 2;
size = rows_p2*cols_p2;
R    = CreateDelOperators(rows_p2, cols_p2);%
dd   = -full(sum(R, 2));%-full(sum(R, 2));%-4*ones(rows * cols, 1);
D    = spdiags(dd, 0, size, size);
invD = spdiags(1./dd, 0, size, size);
[Dx, Dy] = CreateDerivativeOperators(rows_p2, cols_p2);

p = padarray(p, [1 1]);
q = padarray(q, [1 1]);
px = Dx*p(:);
qy = Dy*q(:);
% px = padarray(px, [1 1]);
% qy = padarray(qy, [1 1]);

A = R + D;
b = px + qy;

x0 = zeros(size, 1);
k_max = 5000;
tol = 1e-6;

x = x0;
k = 1;
while( norm(A*x-b) > tol && k <= k_max )
    x = invD*(b-R*x);
    k = k + 1;
end

Z = reshape(x, [rows_p2 cols_p2]);

figure;
colormap gray;
surf(Z,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
axis tight
shading interp
view(110,45)
camlight left;
axis off
title('Estimated Depth image')

figure;
colormap gray;
surf(mozart,'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud')
axis tight
shading interp
view(110,45)
camlight left;
axis off
title('Ground Truth')

