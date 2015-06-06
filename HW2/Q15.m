close all; clear;
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
epsilon = 1e-12;
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
    'Display','final-detailed', ...
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
svm_opts = statset('MaxIter',200000);
SVMStruct = svmtrain(X,y,'METHOD','SMO','options',svm_opts,'BOXCONSTRAINT',ub,'AUTOSCALE',false);
y_est_ref = svmclassify(SVMStruct,X');
accuracy_train_ref = sum(y_est_ref==y)/length(y)*100;

%% Solve with augmented lagrangian - not finished
ones_N = ones(N, 1);
lambda0 = 0.5*ones_N*C;
mu0 = 10;
eta0 = 1;
beta = 10.01;
ytag_y = y'*y;
MaxIterAug = 1000;
%init
mu = mu0;
eta = eta0;
alpha = lambda0;

%iterate
CostTot = zeros(1, MaxIterAug);
CostMin = zeros(1, MaxIterAug);
EqCost = zeros(1, MaxIterAug);
options = optimset('Algorithm', 'trust-region-reflective',...
    'Display','final-detailed', ...
    'LargeScale', 'on', ...
    'TolFun',tol, ...
    'MaxIter', maxIter);
for k = 1:MaxIterAug,
    Htild = He + mu*ytag_y;%For augmented form
    b     = -ones(N, 1) - eta*y;%For augmented form
    f = @(x)(0.5* x' * Htild * x + mu/2 * (y'* x)^2 - ones_N' * x - eta * x' * y);
    gradf = @(x)(Htild* x + mu * y * (y'* x) - ones_N - eta * y) ;
%     [lambda, CostTot(k)] = ProjectedNewton(H, b, lb, ub, lambda, 2000, tol, tolkkt);
%     [lambda, CostTot(k)] = ProjectedNewton_v2(H, f, gradf, b, lb, ub, lambda, 2000, tol, tolkkt);

    [alpha, CostTot(k)] = GradientProjection(f, gradf, lb, ub, alpha, maxIter, tol);
%     [alpha, CostTot(k)] = quadprog(Htild, b, ...
%                      [], [], ...
%                      [], [], ...%equality cond
%                      lb, ub, ...%box constraints
%                      alpha, ...%starting point
%                      options) ;%instead of projected Newton
    EqCost(k) = y'*alpha;
    CostMin(k) = f(alpha);
    %update mu,eta
    eta = eta - mu*y'*alpha;
    mu  = mu*beta;

    disp('Augmented Lagragian iteration')
    disp(k)
    disp('Distance from ref:')
    norm(alpha - lambda)
end
norm(alpha - lambda)
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

