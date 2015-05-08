function [x, cost]=CG_LS1(x0, A, b, tol,MaxIt)
%CG_LS - x_min = argmin 0.5*x'*Q*x + b'*x
% x_min=cg_ls(zeros(25,1),A,y, 1e-13,59);
%
x=x0;

%When minimizing 0.5x'Qx +b'x ,as in CG lecture, pg4 section3, this is why there is factor 2 on Q:
g =  2*(A')*( (A*x) - b); 

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
    g_new = (A')*(A*x_new -  b);
    
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
    
end


disp(['CG converged at ' num2str(k) ' iterations']);


end

