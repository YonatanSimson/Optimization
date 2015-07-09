%% Q18 - Surface normals

load Images.mat
load LightSources.mat
load mozart.mat

Images = double(Images);
L = double(LightSources);
%N - number of pictures
[rows, cols, N] = size(Images);
Image_n = zeros(rows, cols, 3);

Lpinv = pinv(L);
for n = 1:cols,
    for m = 1:rows,
        Imn = squeeze(Images(m, n, :));
        Image_n(m, n, :) = Lpinv*Imn;
    end
end

%find p, q, N = (-p, -q, 1)
p = -Image_n(:, :, 1)./((Image_n(:, :, 3))  + eps);
q = -Image_n(:, :, 2)./((Image_n(:, :, 3))  + eps);


%% Q19 - Jacobi method
rows_p2 = rows + 2; 
cols_p2 = cols + 2;
size = rows_p2*cols_p2;
R    = CreateDelOperators(rows_p2, cols_p2);%
dd   = -full(sum(R, 2));
D    = spdiags(dd, 0, size, size);
invD = spdiags(1./dd, 0, size, size);
[Dx, Dy] = CreateDerivativeOperators(rows_p2, cols_p2);

p = padarray(p, [1 1]);
q = padarray(q, [1 1]);
px = Dx*p(:);
qy = Dy*q(:);

A = R + D;
b = px + qy;

x0 = zeros(size, 1);
k_max = 50000;
tol = 1e-1;

x = x0;
k = 1;
% x = A\b;
while( norm(A*x-b) > tol && k <= k_max )
    x = invD*(b-R*x);
    k = k + 1;
end

Z = reshape(x, [rows_p2 cols_p2]);

figure;
colormap gray;
surf(Z,'FaceColor','interp',...
   'EdgeColor','none')
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

