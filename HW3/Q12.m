%% Brachistochrone problem q.12
n    = 500;
g    = 9.8;
[X, Y] = meshgrid(1:n, 1:n);
n_y = 1./sqrt(2 * g * Y);


%% Q12

x0 = [1,   1];%xy
x1 = [300, 300];%xy
S0 = 0;
S  = runFSM(nref,x0,S0);

close all;
figure(1);
imagesc(S);axis square;

hold on; plot(x0(1),x0(2),'sb');
hold on; plot(x1(1),x1(2),'sk');

%% Q13
[Gx, Gy] = gradient(S);
alpha_k = 0.9;

%iterate
xOld = inf(length(x0'), 1);
x = x1';%starting from end point working back to initial position
path = x1';
maxIter = 10000;
tol = 1e-6;
f = @(x)(interp2(X, Y, S, x(1), x(2)));
for k = 1:maxIter,
    dx = -interp2(X, Y, Gx, x(1), x(2), 'linear', 0);
    dy = -interp2(X, Y, Gy, x(1), x(2), 'linear', 0);
    d = [dx; dy];
    
    if (norm(x-xOld)<tol*norm(xOld) || norm(x-x0') < 1 )
        break;
    end
    %a_k = arg min(f + a*d)
    %alpha_k = GoldenSectionLineSearch(@(t)f(x+t*d), 0, alpha_k, maxIter, tol);
    %constant step size
    alpha_k = 0.1 / norm(d);
    %update
    xOld = x;
    x = x + alpha_k*d;
    path = [path x];
    
    hold on; plot([xOld(1); x(1)], [xOld(2); x(2)], '-b', 'LineWidth', 2);

end

%% Q14 - analytic solution
%x = 0.5*k^2*(t - sin(t))
%y = 0.5*k^2*(1 - cos(t))
syms k t0 t1
%find k and t1
sol = solve(0.5*k^2*(t1 - sin(t1)) == x1(1), 0.5*k^2*(1 - cos(t1)) == x1(2));
kk = double(sol.k);

t = linspace(0, double(sol.t1), 1000);
x = 0.5*kk^2*(t - sin(t));
y = 0.5*kk^2*(1 - cos(t));

hold on;
plot(x, y, 'r');



