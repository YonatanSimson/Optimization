function Q6and8(h)

%% define LPF parameters
L = 500;
DeltaP = 0.1;
DeltaS = 0.001;

%% Equiriple design of LPF-using LP
wp = 0.26*pi;
ws = 0.34*pi;
wc = 0.30*pi;
w = (0:L)*pi/L;
%weighting matrix
S = diag(1/DeltaP*(w<=wp) + 1/DeltaS*(w>=ws));
d = (w<=wc)';%desired LPF frequency responce

N = length(h);

T1 = fliplr(tril(ones(N)));
s = T1*h;

figure;
stem(0:length(h)-1, s);
title('Step responce')
ylabel('s[n]');
xlabel('n');

%% Q8 - Time constraints
DeltaT = 0.05;


for N1 = 20:-1:1;
    T2 = aToh(N);
    T = T1(1:N1+1, :)*T2;


    [delta, h_new] = EquirippleDesign(S, d, T);
    displayResults(h, 2)
    disp(['delta: ' num2str(delta)])
    disp(['N1: ' num2str(N1)])

    s_new = T1*h_new;

    figure(2);
    stem(0:length(h_new)-1, s, 'b');
    hold on;
    stem(0:length(h_new)-1, s_new, 'g');
    hold off;
    title('Step responce - time constraints')
    ylabel('s[n]');
    xlabel('n');
    legend('Original step responce', 'Step responce - with constraints')
    if ( delta < 1)
        break;
    end
end

figure;
stem(0:length(h_new)-1, abs(s_new-s), 'b');
title('Step responce difference |s_2[n]-s_1[n]|')
ylabel('s[n]');
xlabel('n');


displayResultsComp(h, h_new, 2);

%% Nested Functions
    function [del, h] = EquirippleDesign(S, H_d, T)
        M = (N-1)/2;
        C = [ones(L+1, 1) cos(w'*(1:M))];
        
        A = [ S*C -ones(L+1,1);
             -S*C -ones(L+1,1);
              T   zeros(N1+1,1);
             -T   zeros(N1+1,1);];
        v = 2*DeltaT*ones(N1+1, 1);
        b = [S*H_d; -S*H_d; v; v];
        f = [zeros(M+1,1); 1];
        
        % solve linear program
        options = optimoptions('linprog', 'Algorithm', 'simplex');
        x = linprog(f, A, b, [], [], [], [], [], options);

        del = x(end);
        a = x(1:M+1);
%         h = [0.5*a(M+1:-1:2); a(1); 0.5*a(2:M+1)];
        T2 = aToh(N);
        h = 0.5*T2*a;
        
    end

end
