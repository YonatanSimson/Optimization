function x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon)

m = size(A, 1);
n = size(A, 2);

%check to see that the intial point is strictly feasible
if ( sum(A*x0 < b) < m )
    x = nan(n, 1);
    return; 
end

t = t0;
x = x0;
maxIter = 10000;
tol = 1e-12;
for k = 1:1000,
    if (m/t < epsilon)
        break;
    end
%     f = @(x)(t*c'*x - sum(log(b - A * x)));
%     grad_f = @(x)(t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2));
%     hessian_f = @(x)(A' * diag((b - A * x).^2)*A);
%     [x, ~] = NewtonMethod_v2(f, grad_f, hessian_f, x, maxIter, tol);
    
    [x, ~] = NewtonMethod(A, b, c, t, x, maxIter, tol);
    disp(['Cost is now: ' num2str(c'*x)]);
    if ( sum(A*x < b) < m )
        error('Point returned from newton step is not strictly feasible'); 
    end

    t = t*mu;
end


