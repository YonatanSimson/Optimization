function [x, Cost] = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon)

m = size(A, 1);

%check to see that the intial point is strictly feasible
if ( sum(A*x0 < b) < m )
    error('Initial point is infeasible')
end

t = t0;
x = x0;
maxIter = 1000;
tol = 1e-6;
OldCost = inf;

for k = 1:100000,
    if (m/t < epsilon)
        disp('m/t < epsilon')
        break;
    end
%     options = optimoptions(@fminunc,'GradObj','on','Hessian','on', ...
%         'MaxIter', maxIter, 'tolX', tol, 'MaxFunEvals', maxIter);%, 'DerivativeCheck', 'on');
%     [x_ref, Cost] = fminunc(@myfun, x, options);
    [x, ~] = NewtonMethod(A, b, c, t, x, maxIter, tol);
%     disp(['Diff between reference my func: ' num2str(norm(x-x_ref))])
    Cost = c'*x;
    disp(['LP Cost is now: ' num2str(Cost)]);
    if ( sum(A*x < b) < m )
        error('Point returned from newton step is not strictly feasible'); 
    end
    if (abs(Cost-OldCost) < tol*Cost)
        disp('Cost not changing')
        break;
    end
    OldCost = c'*x;
    t = t*mu;
end

Cost = c'*x;

%% Nested function - for fminunc and verification
    function [f_x, g, H] = myfun(x)
        f_x = t*c'*x - sum(log(b - A * x));    % Cost function
        g = t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2);
        H = A' * diag(1./((b - A * x).^2))*A;
    end

end
