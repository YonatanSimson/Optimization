clear; close all;

load('X1');
load('X2');
load('X3');

Grad = CalcGrad(X1);

figure(1);
subplot(121);
imagesc(X1); axis image; colormap gray
subplot(122);
imagesc(Grad); axis image; colormap gray
