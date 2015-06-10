function y_est = svm_est(X, w, w0)

y_est = sign(w'*X + w0)';
y_est(y_est==0) = 1;
