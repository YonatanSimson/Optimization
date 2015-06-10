function [ w, w0 ] = GetW( X, y, lambda, tolkkt, C )
%GetW Calculate w, w_0
A = bsxfun(@times, X, y');

%Find Inactive set, 0<lambda<C
bndind = find(lambda > tolkkt * C & lambda < (1 - tolkkt) * C) ;
%Find w
w = A(:, bndind)*lambda(bndind);
%find w0
w0 = mean(y(bndind)'-w'*X(:, bndind));

end

