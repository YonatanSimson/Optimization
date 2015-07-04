%
% [Gx, Gy] = gradient(X);
% Gx = reshape(Dx*X(:), size(X))
% Gy = reshape(Dy*X(:), size(X))

function [Del] = CreateDelOperators(X_rows, X_cols)

%% formulate the Dx, Dy matrices
%Dy ~ diff(X, 1, 1)
dy = CalcBaseDerivative(X_cols);
Dy = kron(eye(X_rows), dy);

%Dx ~ diff(X, 1, 2)
dx = CalcBaseDerivative(X_rows);
Dx = kron(dx, eye(X_cols));

Del = Dx + Dy;
end


function d_base = CalcBaseDerivative(n)
e = ones(n,1);
d_base = spdiags([e e], [-1 1], n, n);
end