%Q13
clear; close all;
t =  [-4 -3 -2 -1 0 1 2 3 4]';
f1 = [0 0 0 0 0.5 1 1 1 1]';
f2 = [0 0.0025 0.0180 0.1192 0.5000 0.8808 0.9820 0.9975 1]';

Dx = CreateDerivativeOperators1D(length(t));

figure(1);

subplot(221);
plot(t, abs(Dx*f1));
xlabel('x');
title('||Dx*f1(x)||_1');

subplot(222);
plot(t, abs(Dx*f2));
xlabel('x');
title('||Dx*f2(x)||_1');

subplot(223);
plot(t, (Dx*f1).^2);
xlabel('x');
title('||Dx*f1(x)||_2');

subplot(224);
plot(t, (Dx*f2).^2);
xlabel('x');
title('||Dx*f2(x)||_2');
