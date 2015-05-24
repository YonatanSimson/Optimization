%Source: http://en.wikipedia.org/wiki/Conjugate_gradient_method
%Solution to x = arg min_x{f(x)}
%f(x) = 0.5*x'*A*x-b'x 
function [x] = ConjGrad(A,b,x, maxIter, tol)
    r=b-A*x;
    p=r;
    rsold=r'*r;
 
    for i=1:maxIter,
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        if sqrt(rsnew)<tol
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
end
