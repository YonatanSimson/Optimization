%% Brachistochrone problem q.10
n    = 500;
g    = 9.8;
Y    = (1:n)';
n_y = repmat(1./sqrt(2 * g * Y), [1 n]);

%cycoloid equations
%y_x =  sin(theta)./(1-cos(theta));

figure;
imagesc(nref)

x0 = [1,1];
x1 = [250, 250];
S0 = 0;
S  = runFSM(nref,x0,S0);

close all;
imagesc(S);axis square;

hold on; plot(x0(1),x0(2),'pb');
hold on; plot(x1(1),x1(2),'pk');


syms x y
S = solve(x^2*y^2 - 2*x - 1,x^2 - y^2 - 1)