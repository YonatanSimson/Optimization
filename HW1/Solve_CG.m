
%clear all
%ToyExample_Init_Script;

% x_min=ToyExample_Solve_CG_b(zeros(25,1),y,A , l, L, 1e-13,100,x_ClosedForm);
% x_min=ToyExample_Solve_CG_b(zeros(25,1),y,A , l, L, 1e-13,100,zeros(25,1));
function [x_min]=Solve_CG(x0,y,A,lambda,L,delta,MaxIt)

bb=-2*(A') *  y;

x=x0;

%When minimizing 0.5x'Qx +b'x ,as in CG lecture, pg4 section3, this is why there is factor 2 on Q:
g= 2*(A')*(A*x)+2*lambda*(L')*(L*x) + bb;  %2*((A')*A+l*(L')*L)*x + bb;   %g=2*Q*x +bb;   %g=Calc_grad(x);

d=-g;

k=0;
while (1)
    
    %1
    alpha = 0.5*-g'*d/ (  (d'*A')*(A*d)+lambda*(d'*L')*(L*d)   )  ;%-g'*d/(d'*2*((A')*A+l*(L')*L)*d);
    
    %2
    x1 = x +alpha*d;
    
    %3
    g1= 2*(A')*(A*x1)+2*lambda*(L')*(L*x1) + bb;  %2*((A')*A+l*(L')*L)*x1 + bb; %g1=Calc_grad(x1);
    
    %4
    beta = g1'*g1/(g'*g);
    
    %5
    d=-g1+beta*d;
    
    if   (k>MaxIt) || (norm(x - x1)<delta)
        break;
    else
        x=x1;
        g=g1;
        k=k+1;
    end
    
end

%Solution Quality:
disp('Number of iterations')
k

x_min=x1;

end