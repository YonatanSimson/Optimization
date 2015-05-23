function alpha = ArmijoRule(f, x_k, f_k, gradf_k, d_k, sigma, beta, alpha_0, lb, ub, flag)
%INPUT:
% f - handle to function f()
% x_k - x at current iteration
% f_k - f(x_k)
% gradk_k - grad(f(x_k))
% d_k - usually equals -grad(f(x_k))
% sigma - a bigger sigma make the search more conservative
% beta - smaller beta makes the line search converge faster at the expense
% of accuracy
% lb - lower bound of box constraints
% ub - upper bound of box constraints
% flag - 0: bactrack on straight line, 1:Project on convex set B
%
%OUTPUT:
% alpha_k - alpha_k ~ arg min(f(x_k + alpa*d_k))
%
%Armijo rule
% Source: http://en.wikipedia.org/wiki/Wolfe_conditions
% f(x_k+alpha*d_k)-f(x_k)<= sigma*grad(f(x_k))'*(x_tild - x_k)
% where x_tild = x_k+alpha*d_k
% Alternatively:
% f(x_k+alpha*d_k)-f(x_k)<= sigma*alpha*grad(f(x_k))'*d_k

%Init
alpha = alpha_0;
x = x_k + alpha*d_k;
%Iterate
while (f(x) - f_k > sigma*gradf_k'*(x - x_k))
    alpha = beta*alpha;%shrink alpha
    %backtrack
    if ( flag == 0 )
        % straight line
        x = x_k + alpha*d_k;
    else
        % projection on convex set B: x = Proj_B(x+alpha*d_k)
        x = Proj_B(x_k + alpha*d_k, lb, ub);%not implemented yet
    end
end







