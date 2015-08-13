function alpha = ArmijoRule_v2(f, x_k, f_k, gradf_k, d_k, sigma, beta, alpha_0)
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
%Armijo rule
% Source: http://en.wikipedia.org/wiki/Wolfe_conditions
% f(x_k+alpha*d_k)-f(x_k)<= sigma*grad(f(x_k))'*(x_tild - x_k)
% where x_tild = x_k+alpha*d_k
% Alternatively - when there is no projection:
% f(x_k+alpha*d_k)-f(x_k)<= sigma*alpha*grad(f(x_k))'*d_k

%Armijo-Goldstein rule
%https://en.wikipedia.org/wiki/Backtracking_line_search

%Init
alpha = alpha_0;
x     = x_k + alpha*d_k;
t     = sigma*d_k'*gradf_k;
k     = 1;
%Iterate
while ( f(x) - f_k > alpha*t )
    alpha = beta*alpha;%shrink alpha
    % straight line
    x = x_k + alpha*d_k;
    %%TODO check if is still interior point

    k = k + 1;
    if ( k >= 100 )
        disp('warning: Armijo rule not converging, k> 100');
    end
end







