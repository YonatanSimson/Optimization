function [x_min]=CG_LS(x0,A,y,tol,MaxIt)
%CG_LS Summary of this function goes here

%Regular run:
% x_min = argmin ||Ax-y||^2
% x_min=cg_ls(zeros(25,1),A,y, 1e-13,59);

x=x0;

%When minimizing 0.5x'Qx +b'x ,as in CG lecture, pg4 section3, this is why there is factor 2 on Q:
g =  2*(A')*( (A*x) -y); 

d=-g;

k=0;
while (1)
    
    %1
    alpha =  0.5*-g'*d/ (  (d'*A')*(A*d)   )  ;
    
    %2
    x1 = x +alpha*d;
        
    %3
    g1 = 2*(A')*( (A*x1) -  y);
    
    %4
    beta = g1'*g1/(g'*g);
    
    %5
    d=-g1+beta*d;
    
    if   (k>MaxIt) || (norm(x - x1)<tol)
        break;
    else
        x=x1;
        g=g1;
        k=k+1;
    end
    
end


x_min=x1;
disp(['CG converged at ' num2str(k) ' iterations']);


end

