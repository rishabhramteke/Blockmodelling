function [val] = jsDivergence(vPmf1, vPmf2)
%
% Jenson-shannon divergence.
%
    vPmfAvg = (vPmf1 + vPmf2) / 2;
    val = 0.5 * klDivergence(vPmf1, vPmfAvg) + 0.5 * klDivergence(vPmf2, vPmfAvg);

end % end of function