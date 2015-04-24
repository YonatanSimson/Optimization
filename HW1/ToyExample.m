load('Y')
[i,j,s] = find(Y);
y = s;

%% Calc A matrix
BoundingBox = [1;5;1;5];%[start_x, end_x, start_y, end_y];
X_rows = (BoundingBox(4) - BoundingBox(3) + 1);
X_cols = (BoundingBox(2) - BoundingBox(1) + 1);
A_cols = X_rows*X_cols;
A_rows = length(Y);
A = zeros(A_rows, A_cols);

%first row (5, 1)
startXY = [2, 5];
direction = -3*pi/4;%up, left
A(1, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%second row (2,2)
startXY = [1, 2];
direction = pi/4;%down, right
A(2, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%Third row (6, 2)
startXY = [4, 1];
direction = pi/2;%going down
A(3, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%4th row (3, 3)
startXY = [1, 4];
direction = 0;%going right
A(4, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%5th row (7, 3)
startXY = [2, 1];
direction = pi/4;%going down, right
A(5, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%6th row (2, 4)
startXY = [1, 2];
direction = 0;%right
A(6, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%7th row (1, 5)
startXY = [1, 1];
direction = 0;%right
A(7, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%8th row (4, 5)
startXY = [1, 5];
direction = -pi/4;%up, right
A(8, :) = calcEmitterToReceiverWeights(startXY, direction, BoundingBox, A_cols, X_rows, X_cols);

%% formulate the L matrix
[Dx, Dy] = CreateDerivativeOperators(X_rows, X_cols);

%% LS solution
A= sparse(A);
lambda = 1e-5;

L = [Dx; Dy];
AA = [A; sqrt(lambda)*L];

B = [y; zeros(size(L, 1), 1)];

Xest = (A'*A + lambda*L'*L)\(A'*y);
Xest = reshape(Xest, [X_rows, X_cols]);

I = 255*exp(-Xest);%Io = 255
figure(1); imagesc(uint8(I)); colormap gray; axis image

%Solution using lsqr
x_ls = lsqr(AA, B, 1e-10, 200);
x_ls = reshape(x_ls, [X_rows, X_cols]);
residual_ls = sum(sum((x_ls - Xest).^2));


%% Amit's solutions
tmp = ((A')*A+lambda*(L')*L)\(A');%inv((A')*A+l*(L')*L) *(A');
x_a = tmp*y;

b_lsqr=(A') *  y;
A_lsqr=((A')*A+lambda*(L')*L);
R = chol(A_lsqr);%A_lsqr=R'*R; Cholesky
x_lsqr=lsqr(A_lsqr,b_lsqr,1e-10,500,R',R);%seems that max iteration num is 25 -> limitted to the size of x

x_pcg = pcg(A_lsqr,b_lsqr,1e-10,500,R',R);
 
%Notes:
norm(x_a-x_lsqr)
norm(x_a-x_pcg)
norm(x_a-x_ls(:))

%% Solve the same problem iteratively
[J, grad] = costFunction(Xest(:), A, L, y, lambda);
% Check gradient
costFunc = @(x) costFunction(x, A, L, y, lambda);
x_0 = 0.5*zeros(size(A, 2), 1);

[~, grad] = costFunc(x_0);
numgrad = computeNumericalGradient(costFunc, x_0);
diffGrad = norm(numgrad-grad)/norm(numgrad+grad);
disp(['difference between pre-calculated and numerical gradient is: ' num2str(max(diffGrad))])

%my iterative solution using steepest descent
numOfIter = 5000;
cost = zeros(1, numOfIter);
alpha = 0.001;
alphaMax = 0.1;
rng(123);
x_0 = 0.5*rand(5*5, 1);
x = x_0;
for k = 1:numOfIter,
    [cost(k), grad] = costFunction(x, A, L, y, lambda);
    if (norm(grad) < 1e-12)
        cost = cost(1:k);
        break;
    end
    %line search for alpha_k
    %minimize f(x_k + a_k*d_k) as a function of a_k
%     alpha_k = BisectionLineSearch(costFunc, 0, alphaMax, x, -grad, 200, 1e-8);
    alpha_k = GoldenSectionLineSearch(@(alpha)costFunc(x-alpha*grad), 0, alpha_k, 200, 1e-8);
    x = x - alpha_k*grad;
end

figure(2); plot(1:length(cost), cost); xlabel('n iteration'); ylabel('cost')
x_gs = reshape(x, [X_rows, X_cols]);

residual_gs = norm(x_ls(:)-x_gs(:));

[x, cost] = GradientDescent(AA, B, x_0, 10000, 1e-6);
figure(3); plot(1:length(cost), cost); xlabel('n iteration'); ylabel('cost')
x_gs2 = reshape(x, [X_rows, X_cols]);
residual_gs2 = sum(sum((x_gs2 - Xest).^2));

 %  Run fminunc to obtain the optimal x
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 200);
x_0 = ones(5*5, 1);
[x_min, cost] = ...
	fminunc(@(t)(costFunction(t, A, L, y, lambda)), x_0, options);
x_min = reshape(x_min, [X_rows, X_cols]);

residual_min = sum(sum((x_min - Xest).^2))


