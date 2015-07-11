% q17
close all
clear all

load('I.mat');

figure;
imagesc(I);
axis image
colormap gray
title('Shaded image')

eps = 1E-10;
F   = sqrt(1./I.^2 - 1 );
Fe  = F + eps * (F == 0);

x0 = [121, 143];
z0 = 0;
z  = runFSM(Fe,x0,z0);

z = -z;
z = (z - min(z(:)))/(max(z(:)) - min(z(:)));

figure; 
surf(z); 
colormap gray
shading interp;
axis('tight');
view(110,45);
axis('off');
camlight
title('|\nablaz|=F(x,y)')


