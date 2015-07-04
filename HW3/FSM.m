function [distMap] = FSM(distMap,G11,G22,G12,M,iter)
%{
Input -
distMap = an initialized distmap with 0 in the source points and "inf" (or some large number ) in all other points
G11,G12,G22 = Metric tensor coefficients
M = matrix of dupdate directions
iter = number of iterations
gpuflag = a flag for the gpu algorithm

Output -
distMap = the distance map
%}


m=size(G11,1);
n=size(G11,2);

for i=1:iter


    %first direction
    distMap=distMap.';
    [distMap]=UpdateFSM(distMap,M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),G11.',G22.',G12.',n,m);
    distMap=distMap.';

    %second direction
    distMap=(flipud(distMap.'));
    [distMap]=UpdateFSM(distMap,M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),flipud(G11.'),flipud(G22.'),flipud(G12.'),n,m);
    distMap=flipud(distMap).';

    %third direction
    [distMap]=UpdateFSM(distMap,M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),G11,G22,G12,m,n);
    
    %fourth direction    
    distMap=flipud(distMap);
    [distMap]=UpdateFSM(distMap,M(1),M(2),M(3),M(4),M(5),M(6),M(7),M(8),flipud(G11),flipud(G22),flipud(G12),m,n);
    distMap=flipud(distMap);

end




