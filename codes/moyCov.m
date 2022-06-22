function [mY, C ,Yc]=moyCov(Y)
 % function for calculating cov and mean vector
 mn = size(Y,2);
 mY  = (sum(Y,2)/mn);
 C = cov(Y');  %sum((Y-mY)*(Y-mY)')/mn
 Yc = minus(Y , mY);
end