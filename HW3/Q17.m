% q17
close all
clear all

Data = load('I.mat');
I    = Data.I ;
clear Data

Data =  load('mozart.mat');
mozart = Data.mozart ;
clear Data

eps = 1E-10;
F       = sqrt(1./I.^2 - 1 );
Fe    = F + eps * (F == 0);

x0 = [128,145];
z0 = 0;
z  = runFSM(Fe,x0,z0);

z11    = z(1,1);
z      = -z + z11;

figure; 
surf(z); 
colormap gray
shading interp;
axis('tight');
view(110,45);
axis('off');
camlight

figure; 
surf(mozart); 
colormap gray
shading interp;
axis('tight');
view(110,45);
axis('off');
camlight

