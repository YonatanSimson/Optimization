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

Xest = inv(A'*A + lambda*L'*L)*A'*y;
Xest = reshape(Xest, [X_rows, X_cols]);

I = 255*exp(-Xest);%Io = 255
figure(1); imagesc(uint8(I)); colormap gray; axis image

%% Solve the same problem iteratively
[J, grad] = costFunction(Xest(:), A, L, y, lambda);
checkGradients(A, L, y, lambda);%numerical check of gradient

%my iterative solution using steepest descent
numOfIter = 400;
cost = zeros(1, numOfIter);
alpha = 0.05;
x_0 = 0.5*ones(5*5, 1);
x = x_0;
for k = 1:numOfIter,
    [cost(k), grad] = costFunction(x, A, L, y, lambda);
    x = x - alpha*grad;
end

figure(2); plot(1:numOfIter, cost); xlabel('n iteration'); ylabel('cost')
x_gs = reshape(x, [X_rows, X_cols]);

residual = sum(sum((x_gs - Xest).^2));
 
 %  Run fminunc to obtain the optimal x
%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);
x_0 = 0.5*ones(5*5, 1);
[x_opt, cost] = ...
	fminunc(@(t)(costFunction(t, A, L, y, lambda)), x_0, options);

x_opt = reshape(x_opt, [X_rows, X_cols]);

Q = (A'*A) + lambda*(L'*L);
b = 2*A'*y;
x_tag = lsqr(Q, b);
x_tag = reshape(x_tag, [X_rows, X_cols]);