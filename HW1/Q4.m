clear; close all;

load('X1');
load('X2');
load('X3');

Grad1 = CalcGrad(X1);
Grad2 = CalcGrad(X2);
Grad3 = CalcGrad(X3);

figure(1);
subplot(121);
imagesc(X1); axis image; colormap gray
title('X1')
subplot(122);
imagesc(Grad1); axis image; colormap gray
title('$ ||\nabla X1||_2 $','interpreter','latex')

figure(2);
subplot(121);
imagesc(X2); axis image; colormap gray
title('X2')
subplot(122);
imagesc(Grad2); axis image; colormap gray
title('$ ||\nabla X2||_2 $','interpreter','latex')

figure(1);
subplot(121);
imagesc(X3); axis image; colormap gray
title('X3')
subplot(122);
imagesc(Grad3); axis image; colormap gray
title('$ ||\nabla X3||_2 $','interpreter','latex')
