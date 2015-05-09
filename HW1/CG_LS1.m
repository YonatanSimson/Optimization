function [x, cost]=CG_LS1(x0, A, b, tol,MaxIt)
%CG_LS - x_min = argmin 0.5*x'*Q*x + b'*x
% x_min=cg_ls1(zeros(25,1),A,y, 1e-13,59);
%
x=x0;

%minimizing 0.5x'A'Ax-b'x ,as in CG lecture, pg4 section3:
s=A*x-b;
g =  (A')*s; %(A')*( A*x - b); 

d=-g;

cost = zeros(1, MaxIt);
for k = 1:MaxIt,
    %0
    r = (A*d);
    
    %1
    alpha =  -g'*d / ( r'*r );
    
    %2
    x_new = x + alpha*d;
        
    %3
    s1= s+alpha*r;     %s1= A*x1-y =  s+alpha*A*d = s+alpha*r
    g_new = (A')*s1;   %(A')*(A*x_new -  b);
    
    %4
    beta = g_new'*g_new/(g'*g);
    
    %5
    d=-g_new + beta*d;
    
    if norm(x - x_new)<tol || norm(g_new)<tol,
        cost = cost(1:k);
        break;
    end
    
    g = g_new;
    x = x_new;
    s=s1;
    
end


disp(['CG converged at ' num2str(k) ' iterations']);


end

