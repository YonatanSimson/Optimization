  
 %ToyExample_Solve_NewtonStep([-500;-500])
%ToyExample_Solve_NewtonStep(zeros(25,1))
function [temp]=ToyExample_Solve_NewtonStep(x0)%x0 = [-1000;-1000];
global Q;
global bb;
global Const;

if 0
    Q=[1 1;
        1 2;];
    bb= [75; 50];
    Const = 100;
    %b1=bb(1);b2=bb(1); a=Q(1,1);d=Q(2,2),c=Q(1,2)
    %x_min=0.5*(d*b1-c*b2)/(c*c-a*d)
    %y_min=0.5*(c*b1-d*b2)/(c*c-a*d)
    
    armjio = 1000;
    delta=1e-5;
else
    ToyExample_Init_Script;
    Q=((A')*A+l*(L')*L);
    bb=(A') *  y;
    Const = 0;
    
    armjio = 1;
    delta=1e-12;
end


for iii=1:1000

    %1. Find direction:
    d = Calc_d(x0);
    d=d/norm(d);
    
    %2. + 3.  Find Step size & update x1
    b=x0 + armjio*d;
    x1 = GSLS(x0,b,delta);
    if 0    
        Fx(x1);
        alpha = norm(x0-x1); %Step size
    end
    
    if (norm (x0-x1) < 10*delta)
        'got to break'
        iii
        break;
    end
    
    if (iii<1000)
        x0=x1; 
    end
end

%Solution Quality:
'My Newton + Golden Quality:'
norm(x_ClosedForm-x1)
end

function [fx] = Fx(x)
global Q;
global bb;
global Const;

fx  = x'*Q*x +bb'*x+ Const;%f = x'Qx  %        (x,y)Q(x,y)'=x^2 + y^2;
end

%calc Newton Step
%x,d are column vector in Rn, n=2 or 25 in this ToyExample
function [d] = Calc_d(x)
global Q;
global bb;

H= 2*Q;
g= 2*Q*x + bb;

%d = inv(H)*(-g);
d= H\(-g);

end

% How to determine 'a' and 'b' !!!???:
%-> a is x0 !?
% find b using armjio?
%In case of convex, can find b using a rough check for b, such that f(b) > f(x0) !?
function [x] = GSLS(a,b,delta)
%tau,delta are scalars
%a,b are start & end points of the line search

if (norm(b-a) < delta)
    error('b-a < delta');
end

%prolog

tau = (3-sqrt(5))/2;%~0.382  lecture!     %tau_wikipedia ??? =(-1 + sqrt(5))/2; % ~0.618 Wikipedia!
p=a+tau*(b-a);
q=b-tau*(b-a);

fa= Fx(a);
fb=Fx(b);
fp = Fx(p);
fq=Fx(q);

N = 1000;
%Search Loop:
for (ii=1:N)
    if  (norm(b-a) < delta)
        %'b-a < delta'
        ii;
        
        a;
        b;
        p;
        q;
        
        Fx(p);
        Fx(q);
        x=(p+q)/2;
        break;
    end
    
    if (fp<fq)
        b=q; fb=fq;
        q=p; fq=fp;
        p=a+tau*(b-a);  fp=Fx(p);%d_p =  p*d; fp=Fx(d_p);
    else
        a=p; fa=fp;
        p=q;fp=fq;
        q=b-tau*(b-a); fq =  Fx(q);%d_q = q*d; fq =  Fx(d_q);
    end
end %for

end

