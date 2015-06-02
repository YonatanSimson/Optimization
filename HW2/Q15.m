%% INPUTS
load('xForTraining.mat')
load('labelsForTraining.mat')
load('coeff.mat')
X = ExtractFeatures(xForTraining, coeff);
wImg  =  sqrt(size(xForTraining, 1));
y = labelsForTraining;

y(y==0) = -1;
y(y==9) = +1;

%% PARAMTERS
C = 1;
epsilon = 1e-6;
maxIter = 1000;%For Projected Newton
tol = 1e-8;%For Projected Newton
n = size(A, 1);
N = size(A, 2);%number of training samples

%% INIT
lb = zeros(N, 1);
ub = C*ones(N, 1);
b  = -ones(N, 1);
A = bsxfun(@times, X, y');
H = A'*A;
He = H + eye(size(H))*epsilon;

%% SOLVE
%Solve dual form using quadprog
K = X'*X ;
n = numel(y) ;
Y = diag(y) ;
options = optimset('Algorithm', 'interior-point-convex',...
    'Display','on', ...
    'LargeScale', 'on', ...
    'TolFun',tol, ...
    'MaxIter', maxIter);
lambda = quadprog(Y*K*Y, - ones(N,1), ...
                 [], [], ...
                 y', 0, ...
                 zeros(N,1), C * ones(N,1), ...
                 [], options) ;


%Solve with augmented lagrangian
% [x, active, Cost] = ProjectedNewton(H, b, lb, ub, maxIter, tol);

%Find Inactive set, 0<lambda<C
bndind = find(alpha > tol * C & alpha < (1 - tol) * C) ;
%Find w
w = A(:, bndind)*lambda(bndind);
%find w0
w0 = mean(y(bndind)-w'*X(:, bndind));

%try on test set


