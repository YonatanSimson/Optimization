function Q5

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


for N = 5:2:101,
    [delta, h] = EquirippleDesign(S, d);
    displayResults(h, 2)
    disp(['N: ' num2str(N)])
    disp(['delta: ' num2str(delta)])
    if ( delta < 1 )
        break;
    end
end

%% Nested Functions
    function [del, h] = EquirippleDesign(S, H_d)
        M = (N-1)/2;
        C = [ones(L+1, 1) 2*cos(w'*[1:M])];
        
        A = [ S*C -ones(L+1,1);
             -S*C -ones(L+1,1);];
        b = [S*H_d; -S*H_d];
        f = [zeros(M+1,1); 1];
        
        % solve linear program
        options = optimoptions('linprog', 'Algorithm', 'simplex');
        x = linprog(f, A, b, [], [], [], [], [], options);

        del = x(end);
        a = x(1:M+1);
        h = [a(M+1:-1:2); a(1); a(2:M+1)];
    end

    function displayResults(hh, fsamp)
        figure(1);
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






