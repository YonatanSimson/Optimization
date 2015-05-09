%IRLS

function [x, cost] = IRLS_grad(A, L, y, x_0, epsilon, lambda, maxIter, tol)
m = size(L, 1);
n = size(L, 2);

x = x_0;
W = speye(m);
cost = zeros(maxIter, 1);
for k = 1:maxIter,
    e = (A*x - y);
    cost(k) = 0.5*(e'*e) + lambda*sum(abs(L*x));
    if (k > 1 && cost(k) > cost(k-1))
        error('IRLS diverges')
    end
    grad    = (A'*(A*x - y) + lambda*L'*W*(L*x));
    if (norm(grad) < tol)
        cost = cost(1:k);
        break;
    end
    
    alpha_k = GoldenSectionLineSearch(@(t)f(x-t*grad), 0, 0.99, maxIter, tol);
    x = x - alpha_k*grad;
    gamma = L*x;
    gamma(abs(gamma)<epsilon) = epsilon;
    W = sparse(1:m,1:m,abs(1./gamma));
end


%% Nested function - for cost evaluation
    function val = f(x)
        e = (A*x - y);
        val = 0.5*(e'*e) + lambda*sum(abs(L*x));
    end

end



