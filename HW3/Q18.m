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
Zx = -Image_n(:, :, 1)./(Image_n(:, :, 3) + eps);
Zy = -Image_n(:, :, 2)./(Image_n(:, :, 3) + eps);


%% Q19 - Jacobi method
R = CreateDelOperators(rows, cols);%
D = -4*spdiags(ones(rows * cols, 1), 0, rows*cols, rows*cols);
invD = -0.25*spdiags(ones(rows * cols, 1), 0, rows*cols, rows*cols);
[Dx, Dy] = CreateDerivativeOperators(rows, cols);

px = Dx*Zx(:);
qy = Dy*Zy(:);
A = R + D;
b = px + qy; 

x0 = zeros(rows*cols, 1);
k_max = 5000;
tol = 1e-6;

x = x0;
k = 1;
while( norm(A*x-b) > tol && k <= k_max )
    x = invD*(b-R*x);
    k = k + 1;
end

Z = reshape(x, [rows cols]);

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

