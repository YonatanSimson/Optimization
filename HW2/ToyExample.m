%Toy Example for Armijo rule


f = @(x)sum((x-[-1;0.5]).^2);
f([0;0])

gradf = @(x)(2*(x-[-1;0.5]));

%params
maxIter = 100;
tol = 1e-8;
sigma = 0.3;
beta  = 0.6;
flag = 1;
%Init
x = [0; 0];
alpha_k = 0.9;
lb = [0; 0];
ub = [1; 1];
%iterate
xOld = [inf; inf];
for k = 1:maxIter,
    d = -gradf(x);
    if (norm(x-xOld)<tol)
        disp(['Converged at iteration ' num2str(k)]);
        break;
    end
    %a_k = arg min(f + a*d)
    %alpha = GoldenSectionLineSearch(@(t)f(x+t*d), 0, alpha_k, maxIter, tol);
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, lb, ub, flag);
    %update
    xOld = x;
    x = Proj_B(x + alpha_k*d, lb, ub);
end


