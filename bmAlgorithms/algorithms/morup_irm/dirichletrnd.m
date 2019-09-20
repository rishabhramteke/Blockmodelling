%-------------------------------------
function y = dirichletrnd(a)
% Sample from Dirichlet distribution.
%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986

% This is essentially a generalization of the method for Beta rv's.
% Theorem 4.1, p.594

y = gamrnd(a,1);
y=y/sum(y);