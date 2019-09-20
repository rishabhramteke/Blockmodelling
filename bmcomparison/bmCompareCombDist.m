function [combDist] = bmCompareCombDist(posDist, reconDist, sCombDist, param)
%
% Compute various combing distances.
%

switch sCombDist
    case 'linear'
        % linear weighted sum
        % param here is alpha
        combDist = (1 - param) * posDist + param * reconDist;
    case 'minkowski'
        % minkowski
        % param here is p
        combDist = (posDist^param + reconDist^param)^(1/param);
    otherwise
        error('Unknown combine distance specified')
end

end