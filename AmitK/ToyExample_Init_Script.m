
A = zeros(8,25);
y=zeros(8,1);
X = zeros(5,5);

%Assumig that the length of a cell side is 1 (and all sides equale).   -> y=A*X(:);%y=A*vec(X);

% y(1,1)= 5->1 (from source 5 to reciver 1)
% y(2,1)= 2->2
% y(3,1)= 6->2
% y(4,1)= 3->3
% y(5,1)= 7->3
% y(6,1)= 2->4
% y(7,1)= 1->5
% y(8,1)= 4->5

 [Dx, Dy] = CreateDerivativeOperators(5, 5);
 L=[Dx;Dy];
 l=1e-5;%lambda
 d=sqrt(2);
A=...
[ ...
0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;        % 5->1 (from source 5 to reciver 1)
0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0;      % 2->2
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0;         % 6->2
0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0;         % 3->3
0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0 0 0 0 0 d 0;      % 7->3
0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0;          % 2->4
1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0;          % 1->5
0 0 0 0 d 0 0 0 d 0 0 0 d 0 0 0 d 0 0 0 d 0 0 0 0;...   % 4->5
];

% Y =
% 
%    (5,1)       0.2322
%    (2,2)       1.9477
%    (6,2)       1.5065
%    (3,3)       1.5065
%    (7,3)       1.9477
%    (2,4)       1.5065
%    (1,5)       0.3361
%    (4,5)       2.5065
y=[  0.2322   1.9477      1.5065      1.5065      1.9477    1.5065  0.3361    2.5065]';

%Normal equations :  ((A')*A+l*(L')*L)x - (A')*y =0;

%Closed form solution:
tmp=  ((A')*A+l*(L')*L)\(A');%inv((A')*A+l*(L')*L) *(A');
x_ClosedForm=tmp*y;

XX = reshape(x_ClosedForm,5,5);

 %Lsqr solution:
 b_wave=(A') *  y;
 A_wave=((A')*A+l*(L')*L);
 x_lsqr=lsqr(A_wave,b_wave,1e-30,30);%seems that max iteration num is 25 -> limitted to the size of x

 x_pcg = pcg(A_wave,b_wave,1e-30,100);
 
 %Notes:
%  norm([x-x_lsqr])
% 
% ans =
% 
%     0.3999
% 
% norm([x-x_pcg])
% 
% ans =
% 
%   1.1204e-010
 
%Precoditioning:
R = chol(A_wave);%A_lsqr=R'*R; Cholesky
P = inv(R);
x_lsqr_preCond=lsqr(A_wave,b_wave,1e-30,30,(inv(P))',inv(P));
disp('lsqr precond quality:')
norm(x_ClosedForm-x_lsqr_preCond)

x_pcg_precond = pcg(A_wave,b_wave,1e-30,100,(inv(P))',inv(P));
disp('pcg precond quality:')
norm(x_ClosedForm-x_pcg_precond)
 