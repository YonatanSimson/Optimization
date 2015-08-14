close all; clear;
Flag = 'projGrad'; %projGrad, projNewton, quadprog

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
maxIter = 200000;%For Projected Newton/Gradient Projection
tol = 1e-10;%For Projected Newton/Gradient Projection
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
if ( 1 )
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

[ w_qp, w0_qp ] = GetW( X, y, lambda, tolkkt, C );
y_train_qp = svm_est(X, w_qp, w0_qp);

%accuracy for training set
accuracy_train_qp = sum(y_train_qp==y)/length(y)*100;


%% Matlab SVM reference
svm_opts = statset('MaxIter',200000);
SVMStruct = svmtrain(X,y,'METHOD','SMO','options',svm_opts,'BOXCONSTRAINT',ub,'AUTOSCALE',false);
y_train_ref = svmclassify(SVMStruct,X');
accuracy_train_ref = sum(y_train_ref==y)/length(y)*100;
end
%% Solve with augmented lagrangian
ones_N = ones(N, 1);
lambda0 = 0.5*ones_N*C;
mu0 = 10;
eta0 = 1;
beta = 10.01;
y_y_tag = y*y';
MaxIterAug = 1000;
maxIter    = 50000;
tol        = 1e-6;

%init
mu = mu0;
eta = eta0;
alpha = lambda0;
alphaOld = inf(length(lambda0), 1);

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
    if ( sum((y'*alpha).^2) < tol && ...
            norm(alpha-alphaOld) < tol)
        CostTot = CostTot(1:k-1);
        CostMin = CostMin(1:k-1);
        EqCost = EqCost(1:k-1);
        
        disp('Augmented Lagragian converged at iteration')
        disp(k)
        break;
    end
    Htild = He + mu*y_y_tag;%For augmented form
    b     = -ones(N, 1) - eta*y;%For augmented form
    f = @(x)(0.5* x' * He * x + mu/2 * (y'* x)^2 - ones_N' * x - eta * x' * y);
    gradf = @(x)(He* x + mu * y * (y'* x) - ones_N - eta * y) ;
    if ( strcmp(Flag, 'projGrad') )
        [alpha, CostTot(k)] = GradientProjection(f, gradf, lb, ub, alpha, maxIter, tol);
    elseif( strcmp(Flag, 'projNewton') )
        [alpha, CostTot(k)] = ProjectedNewton(Htild, b, lb, ub, alpha, maxIter, tol, tolkkt);
    else
        [alpha,  CostTot(k)] = quadprog(Htild, b, ...
                         [], [], ...
                         [], [], ...%equality cond
                         lb, ub, ...%box constraints
                         alpha, ...%starting point
                         options) ;
    end
    EqCost(k) = (y'*alpha).^2;
    CostMin(k) = 0.5* alpha' * H * alpha - ones_N' * alpha;
    %update mu,eta
    eta = eta - mu*y'*alpha;
    mu  = mu*beta;

    disp('Augmented Lagragian iteration')
    disp(k)
    disp('Distance from ref:')
%     norm(alpha - lambda)
    alphaOld = alpha;
end
% norm(alpha - lambda)

%train results from augmeted lagrangian
[ w_al, w0_al ] = GetW( X, y, alpha, tolkkt, C );
y_train_al = svm_est(X, w_al, w0_al);
%accuracy for training set
accuracy_train_al = sum(y_train_al==y)/length(y)*100;
disp(['Trainining set accuracy: ' num2str(accuracy_train_al)])

%% Test - Q16
load('xForTest.mat')
load('labelsForTest.mat')
Xtest = ExtractFeatures(xForTest, coeff);
y_test = labelsForTest;

y_test(y_test==0) = -1;
y_test(y_test==9) = +1;


%accuracy for test set - AL
y_test_al = svm_est(Xtest, w_al, w0_al);
accuracy_train_al = sum(y_test_al==y_test)/length(y_test)*100;

%accuracy for test set - Reference
y_test_ref = svmclassify(SVMStruct, Xtest');
accuracy_test_ref = sum(y_test_ref==y_test)/length(y_test)*100;

disp(['Test set accuracy: ' num2str(accuracy_train_al)])

