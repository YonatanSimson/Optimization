function alpha = ArmijoRule(f, A, b, x_k, f_k, gradf_k, d_k, sigma, beta, alpha_0)
%INPUT:
% f - handle to function f()
% x_k - x at current iteration
% f_k - f(x_k)
% gradk_k - grad(f(x_k))
% d_k - can be -grad(f(x_k)) or H^-1*grad(f(x_k))
% sigma - a bigger sigma make the search more conservative
% beta - smaller beta makes the line search converge faster at the expense
% of accuracy
%
%OUTPUT:
% alpha_k - alpha_k ~ arg min(f(x_k + alpa*d_k))
%
%Armijo-Goldstein rule
%https://en.wikipedia.org/wiki/Backtracking_line_search

%Init
m     = size(A, 1); 
alpha = alpha_0;
x     = x_k + alpha*d_k;
t     = -sigma*d_k'*gradf_k;
k     = 1;
%Iterate
while (f(x) - f_k > sigma*gradf_k'*(x - x_k) || sum(A*x-b<-eps) < m )%if point not strictly feasible keep on shrinking
    alpha = beta*alpha;%shrink alpha
    % straight line
    x = x_k + alpha*d_k;
    %%TODO check if is still interior point

    k = k + 1;
    if ( k >= 1000 )
        disp('warning: Armijo rule not converging, k> 100');
    end
end







