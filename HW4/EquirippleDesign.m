function [del, h] = EquirippleDesign(S, d, L, N)
    M = (N-1)/2;
    w = (0:L)*pi/L;

    C = [ones(L+1, 1) cos(w'*(1:M))];

    A = [ S*C -ones(L+1,1);
         -S*C -ones(L+1,1);];
    b = [S*d; -S*d];
    f = [zeros(M+1,1); 1];

    % solve linear program
    options = optimoptions('linprog', 'Algorithm', 'simplex');
    x = linprog(f, A, b, [], [], [], [], [], options);
    del = x(1);
    a = x(2:M+2);
    T2 = aToh(N);
    h = 0.5*T2*a;
end


