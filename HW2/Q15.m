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
maxIter = 200000;%For Projected Newton
tol = 1e-8;%For Projected Newton
tolkkt = 1e-3;
n = size(X, 1);%dimension of w
N = size(X, 2);%number of training samples

%% INIT
lb = zeros(N, 1);
ub = C*ones(N, 1);
b  = -ones(N, 1);
A = bsxfun(@times, X, y');
H = A'*A;
He = H + eye(size(H))*epsilon;

%% SOLVE
%Solve dual form using quadprog
n = numel(y) ;
Y = diag(y) ;
%K = X'*X; H = Y*K*Y; Andrea Vadaldi form
options = optimset('Algorithm', 'interior-point-convex',...
    'Display','on', ...
    'LargeScale', 'on', ...
    'TolFun',tol, ...
    'MaxIter', maxIter);
lambda = quadprog(H, b, ...
                 [], [], ...
                 y', 0, ...%equality cond
                 lb, ub, ...%box constraints
                 [], options) ;
             
%Find Inactive set, 0<lambda<C
bndind = find(lambda > tolkkt * C & lambda < (1 - tolkkt) * C) ;
%Find w
w = A(:, bndind)*lambda(bndind);
%find w0
w0 = mean(y(bndind)'-w'*X(:, bndind));

y_est = sign(w'*X + w0)';
y_est(y_est==0) = 1;
%accuracy for training set
accuracy_train = sum(y_est==y)/length(y);


%% Matlab SVM reference
opts = statset('MaxIter',200000);
SVMStruct = svmtrain(X,y,'METHOD','SMO','options',opts,'BOXCONSTRAINT',ub,'AUTOSCALE',false);
y_est_ref = svmclassify(SVMStruct,X');
accuracy_train_ref = sum(y_est_ref==y)/length(y)*100;

%% Solve with augmented lagrangian - not finished
lambda0 = 0.5*ones(N, 1)*C;
mu0 = 1;
eta0 = 1;
beta = 10;
ytag_y = y'*y;

mu = mu0;
eta = eta0;
lambda = lambda0;

%iterate
Htild = He + mu*ytag_y;%For augmented form
b     = ones(N, 1) + eta*y;%For augmented form
[lambda, active, Cost] = ProjectedNewton(H, b, lb, ub, lambda, maxIter, tol, tolkkt);
%update mu,eta
eta = eta - mu*y'*lambda;
mu  = mu*beta;


%% Test
load('xForTest.mat')
load('labelsForTest.mat')
Xtest = ExtractFeatures(xForTest, coeff);
y_test = labelsForTest;

y_test(y_test==0) = -1;
y_test(y_test==9) = +1;


y_test_est = sign(w'*Xtest + w0)';
y_test_est(y_test_est==0) = 1;


accuracy_test = sum(y_test_est==y_test)/length(y_test)*100;

