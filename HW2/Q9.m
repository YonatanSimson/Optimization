load('xForTraining.mat')
load('labelsForTraining.mat')
load('coeff.mat')
x_ = ExtractFeatures(xForTraining, coeff);
w  =  sqrt(size(xForTraining, 1));
y = labelsForTraining;

y(y==0) = -1;
y(y==9) = +1;

At = w'*x_ + w0 < 1;%set of active samples
    
CurrCost = lambda/2*(w'*w) + 1/N*sum((1-(w'*x_ + w0)).*At);

SubGradientW  = 1/N*x_(:, At)*y(At);
SubGradientw0 = 1/N*sum(y(At));
    
