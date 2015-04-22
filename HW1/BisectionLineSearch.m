function alpha_k = BisectionLineSearch(func, a, b, x, d, maxIter, tol)
% Do a line search for a minimum of func = @(alpha)f(x+alpha*d), where
% alpha is in the interval [a,b]. 
% We assume that func also returns the derivative of the function by x

for k = 1:maxIter,
    t = (a + b)/2;
    gTagAlpha = CalcGradientG(t);
    if ( abs(gTagAlpha) < tol )
        break;
    end
    if (gTagAlpha > 0)
        b = t;
    else
        a = t;
    end
end

alpha_k = t;
%define g'(alpha) = grad(f(x+alpha*d))'*d
    function gTagAlpha = CalcGradientG(alpha)
        [~, gradX] = func(x+alpha*d);
        gTagAlpha  = gradX'*d;
    end
end