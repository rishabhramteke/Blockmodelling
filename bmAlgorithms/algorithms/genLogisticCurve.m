function [vX, vY] = genLogisticCurve(xmin, xmax, numXPts, A, K, B, v, Q, M)
%
% Return the generalised logistic curve within the linspace setings.
%
% A - lower asymtotic value
% k - upper asymptotic value
% B - growth rate
% v - near which asumptote is the maximum growth occuring
% Q - depend on Y(0)
% M - time of maximum growth if Q = v
%

vX = linspace(xmin, xmax, numXPts);


vY = A + (K-A) ./ (1 + Q * exp(-B/v * (vX - M)));

plot(vX, vY);

