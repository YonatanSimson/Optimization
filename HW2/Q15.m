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
bndind = find(lambda > tolkkt * C & lambda < (1 - tolkkt) * C) ;
%Find w
w = A(:, bndind)*lambda(bndind);
%find w0
w0 = mean(y(bndind)'-w'*X(:, bndind));

y_est = sign(w'*X + w0)';
y_est(y_est==0) = 1;
%accuracy for training set
accuracy_train = sum(y_est==y)/length(y);

opts = statset('MaxIter',200000);
SVMStruct = svmtrain(X,y,'METHOD','SMO','options',opts,'BOXCONSTRAINT',ones(N, 1)*C,'AUTOSCALE',false);
y_est_ref = svmclassify(SVMStruct,X');
accuracy_train_ref = sum(y_est_ref==y)/length(y)*100;


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

