function [S] = runFSM(nref,x0,varargin)
%{
Fast Sweeping algorithm
Input-
nref - non-negative index of refraction
x0 - a 2x1 vector of coordinates of the source point [x,y]
varargin : S0 - values of the distance map at the points x0; (default : S0=0);

Output-
S - optical length on a regular grid

copyright 2014 Amit Boyarski
%}
if nargin>2
    S0=varargin{1};
else
    S0=zeros(numel(x0),1);
end

    

m=size(nref,1);
n=size(nref,2);
infNum=10^10;
S=zeros(m,n)+infNum;
idx=sub2ind([m,n],x0(:,2),x0(:,1));
S(idx)=S0;
S = padarray(S,[1 1],infNum,'both');

G11=padarray(nref.^2,[1 1],infNum,'both');
G22=G11;
G12=zeros(size(G11));

M=[ -1 -1 0 -1 0 -1 1 -1];

iter=4;
disp('Please wait...');
[S] = FSM(S,G11,G22,G12,M,iter);
S=S(2:end-1,2:end-1);
disp('Done...')

end