function alpha_k = GoldenSectionLineSearch(func, a, b, maxIter, tol)
% Do a line search for a minimum of func = @(alpha)f(alpha), where
% alpha is in the interval [a,b]. 

tau = (3 - sqrt(5))/2;
% fa = func(a);
% fb = func(d);
p = a + tau*(b-a);
q = b - tau*(b-a);
fp = func(p);
fq = func(q);

for k = 1:maxIter,
    if ( abs(b-a) < tol )
        break;
    end
    if (fp < fq)
        b = q; %fb = fq;
        q = p; fq = fp;
        
        p = a + tau*(b-a);
        fp = func(p);
    else
        a = p; %fa = fp;
        p = q; fp = fq;
        q = b - tau*(b-a);
        fq = func(q);
    end
end

alpha_k = (a+b)/2;

