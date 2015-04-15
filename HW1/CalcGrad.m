function Grad = CalcGrad(X)

[X_rows, X_cols] = size(X);
[Dx, Dy] = CreateDerivativeOperators(X_rows, X_cols);
Gx = reshape(Dx*X(:), size(X));
Gy = reshape(Dy*X(:), size(X));


Grad = sqrt(Gx.^2 + Gy.^2);
