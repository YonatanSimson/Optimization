%
% projection onto active set
%
function px = Proj_B(x, lb, ub)
ndim=length(x);
px=zeros(ndim,1); %#ok<NASGU>
px=max(lb, x);
px=min(ub, px);