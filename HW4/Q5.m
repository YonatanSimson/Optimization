function h=Q5

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
    %Once delta is smaller than 1 we know that we satisfy the design
    %constraints
    if ( delta < 1 )
        break;
    end
end

%% Nested Functions
    function [del, h] = EquirippleDesign(S, H_d)
        M = (N-1)/2;
        C = [ones(L+1, 1) cos(w'*[1:M])];
        
        A = [ S*C -ones(L+1,1);
             -S*C -ones(L+1,1);];
        b = [S*H_d; -S*H_d];
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






