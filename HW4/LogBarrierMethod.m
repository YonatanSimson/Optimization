function x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon)

m = size(A, 1);
n = size(A, 2);

%check to see that the intial point is strictly feasible
if ( sum(A*x0 < b) < m )
    x = nan(n, 1);
    error('Initial point is infeasible')
    return; 
end

t = t0;
x = x0;
maxIter = 10000;
tol = 1e-6;
OldCost = inf;
for k = 1:1000,
    if (m/t < epsilon)
        break;
    end
%     f = @(x)(t*c'*x - sum(log(b - A * x)));
%     grad_f = @(x)(t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2));
%     hessian_f = @(x)(A' * diag(1./((b - A * x).^2))*A);
%     [x, ~] = NewtonMethod_v2(f, grad_f, hessian_f, x0, maxIter, tol);
    options = optimoptions(@fminunc,'GradObj','on','Hessian','on', ...
        'MaxIter', maxIter, 'tolX', tol, 'MaxFunEvals', maxIter);%, 'DerivativeCheck', 'on');
    [x_ref, Cost] = fminunc(@myfun, x, options);
    [x, ~] = NewtonMethod(A, b, c, t, x, maxIter, tol);
    disp(['Diff between reference my func: ' num2str(norm(x-x_ref))])
    disp(['LP Cost is now: ' num2str(c'*x)]);
    if ( sum(A*x < b) < m )
        error('Point returned from newton step is not strictly feasible'); 
    end
    if (abs(Cost-OldCost) < tol*Cost)
        disp('Cost not changing')
        break;
    end

    t = t*mu;
end



    function [f_x, g, H] = myfun(x)
        f_x = t*c'*x - sum(log(b - A * x));    % Cost function
%debug code
%         gg = t*c;
%         HH = zeros(n);
%         for l = 1:m,
%             a_l = A(l, :)';
%             scale = 1/(b(l) - a_l' * x);
%             gg = gg + a_l * scale;
%             HH = HH + a_l * a_l' * scale^2;
%         end
        g = t*c + sum(bsxfun(@rdivide, A', (b - A * x)'), 2);
        H = A' * diag(1./((b - A * x).^2))*A;
    end

end
