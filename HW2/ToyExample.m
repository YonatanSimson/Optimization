%Toy Example for Armijo rule


f = @(x)sum((x-[-1;0.5]).^2);
f([0;0])

gradf = @(x)(2*(x-[-1;0.5]));

%params
maxIter = 100;
tol = 1e-8;
sigma = 0.3;
beta  = 0.6;
flag = 0;
%Init
x = [0; 0];
alpha_k = 1;
%iterate
for k = 1:maxIter,
    d = -gradf(x);
    if (norm(d)<tol)
        disp(['Converged at iteration ' num2str(k)]);
        break;
    end
    %a_k = arg min(f + a*d)
    alpha = GoldenSectionLineSearch(@(t)f(x+t*d), 0, alpha_k, maxIter, tol);
    alpha_k = ArmijoRule(f, x, f(x), gradf(x), d, sigma, beta, alpha_k, flag);
    x = x + alpha_k*d;
end


