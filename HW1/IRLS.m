%IRLS

function [x, cost] = IRLS(A, L, y, x_0, epsilon, lambda, maxIter, tol)
m = size(L, 1);
n = size(L, 2);

x = x_0;
W = speye(m);
cost = zeros(maxIter, 1);

%Solution using lsqr
B = [y; zeros(size(L, 1), 1)];
x_cg = x_0;
for k = 1:maxIter,
    AA = [A; sqrt(lambda)*sqrt(W)*L];
    x_cg = lsqr(AA, B, tol, 500, [], [], x_cg);%replace with CG

    cost(k) = f(x_cg);
    if (k > 1 && cost(k) > cost(k-1))
        error('IRLS diverges')
    end
    
    gamma = L*x_cg;
    gamma(abs(gamma)<epsilon) = epsilon;
    W = sparse(1:m,1:m,abs(1./gamma));
end



%% Nested function - for cost evaluation
    function val = f(x)
        e = (AA*x - B);
        val = 0.5*(e'*e) + lambda*sum(abs(L*x));
    end

end


