function x = LogBarrierMethod(A, b, c, t0, x0, mu, epsilon)

m = size(A, 1);
n = size(A, 2);

%check to see that the intial point in feasible
if ( sum(A*x0 <= b) < m )
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
    [x, ~] = NewtonMethod(A, b, c, t, x, maxIter, tol);
    t = t*mu;
end


