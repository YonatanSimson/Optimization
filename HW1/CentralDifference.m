function Xtag = CentralDifference(X, dim)

[m, n] = size(X);
if ( dim == 1 )%Dy
    Xtag = zeros(m-1, n);
    Xtag(1,:) = X(2, :) - X(1, :);
    Xtag(m-1,:) = X(m, :) - X(m-1, :);
    for k = 2:m-2,
        Xtag(k,:) = 0.5*(X(k+1, :) - X(k-1, :)); 
    end
elseif (dim == 2)%Dx
    Xtag = zeros(m, n-1);
    Xtag(:, 1)   = X(:, 2) - X(:, 1);
    Xtag(:, n-1) = X(:, n) - X(:, n-1);
    for k = 2:n-2,
        Xtag(:, k) = 0.5*(X(k+1, :) - X(k-1, :)); 
    end
else
    error('Illegal value for ''dim''');
end
    