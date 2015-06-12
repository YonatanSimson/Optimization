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
MaxIterAug = 1000;
maxIter    = 50000;
tol        = 1e-6;
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

%% Solve with augmented lagrangian
ones_N = ones(N, 1);
lambda0 = zeros(N, 1);%0.5*ones_N*C;
mu0 = 10;
eta0 = 1;
beta = 10.01;
y_y_tag = y*y';

%init
mu = mu0;
eta = eta0;
lambda = lambda0;
lambdaOld = inf(length(lambda0), 1);

%iterate
CostTot = zeros(1, MaxIterAug);
options = optimset('Algorithm', 'trust-region-reflective',...
    'Display','final-detailed', ...
    'LargeScale', 'on', ...
    'TolFun',tol, ...
    'MaxIter', maxIter);
tic;
for k = 1:MaxIterAug,
    if ( sum((y'*lambda).^2) < tol && ...
            norm(lambda-lambdaOld) < tol*norm(lambda))
        CostTot = CostTot(1:k-1);
        
        break;
    end
    Htild = He + mu*y_y_tag;%For augmented form
    b     = -ones(N, 1) - eta*y;%For augmented form
    f = @(x)(0.5* x' * He * x + mu/2 * (y'* x)^2 - ones_N' * x - eta * x' * y);
    gradf = @(x)(He* x + mu * y * (y'* x) - ones_N - eta * y) ;
    if ( strcmp(Flag, 'projGrad') )
        [lambda, CostTot(k)] = GradientProjection(f, gradf, lb, ub, lambda, maxIter, tol);
    elseif( strcmp(Flag, 'projNewton') )
        [lambda, CostTot(k)] = ProjectedNewton(Htild, b, lb, ub, lambda, maxIter, tol, tolkkt);
    else
        [lambda,  CostTot(k)] = quadprog(Htild, b, ...
                         [], [], ...
                         [], [], ...%equality cond
                         lb, ub, ...%box constraints
                         lambda, ...%starting point
                         options) ;
    end
    %update mu,eta
    eta = eta - mu*y'*lambda;
    mu  = mu*beta;

    lambdaOld = lambda;
    disp(['Augmented Lagragian iteration: ' num2str(k)])
end
toc
%train results from augmeted lagrangian
[ w_al, w0_al ] = GetW( X, y, lambda, tolkkt, C );
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

disp(['Test set accuracy: ' num2str(accuracy_train_al)])

