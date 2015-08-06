function Q5

%% define LPF parameters
N = 101;
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

% [del, hOpt] = EquirippleDesign(S, d, L, N);

% figure;
% plot()
% C = [ones(L+1, 1) cos(w'*(1:M))];
% HOpt = C*hOpt;
% H_d  = 
% displayResults(hOpt, pi)


[delta, h] = EquirippleDesign(S, d);

displayResults(h, pi)
disp(['delta: ' num2str(delta)])

%% Nested Functions
    function [del, h] = EquirippleDesign(S, H_d)
        M = (N-1)/2;
        C = [ones(L+1, 1) 2*cos(w'*[1:M])];
        
        A = [-ones(L+1,1)  S*C;
            -ones(L+1,1) -S*C;];
        b = [S*H_d; -S*H_d];
        f = [1; zeros(M+1,1)];
        
        % solve linear program
        %x = linprog(f,A,b);
        options = optimoptions('linprog', 'Algorithm', 'simplex');
        x = linprog(f, A, b, [], [], [], [], [], options);

        del = x(1);
        a = x(2:M+2);
        h = [a(M+1:-1:2); a(1); a(2:M+1)];
    end

    function displayResults(hh, fsamp)
        figure;
        subplot(2,1,1)
        [H,f] = freqz(hh,1,1024,fsamp);
        semilogy(f,abs(H)), grid on
        title('Frequency responce F\{h_d[n]\}')
        ylabel('|H_d(\omega)|')
        xlabel('Frequency[Hz]')

        subplot(2,1,2)
        stem(0:length(hh)-1, hh);
        title('Impulse responce')
        ylabel('h_d[n]');
        xlabel('n');

    end
end






