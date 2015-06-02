function [x_] = ExtractFeatures(x,coeff)
%{
Extracts 50 principle components of x
Input - x - vector of raw data
        coeff - basis tranformation matrix

Output - x_ - 50 lragest principle components of x
%}
NumofPrincComp=50;
x_=coeff\x;
x_=x_(1:NumofPrincComp,:);

end

