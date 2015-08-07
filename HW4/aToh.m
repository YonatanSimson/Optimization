function T2 = aToh(N)
    M = (N-1)/2;

    T2 = eye(M+1);
    T2 = [flipud(T2(2:end, :)); T2];
    T2(M + 1, 1) = 2;
end