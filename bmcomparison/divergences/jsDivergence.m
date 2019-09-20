function div =  jsDivergence(vPmf1, vPmf2, vWeightDistValues)
%
% Computes the JS divergence between the two pmfs.
%
% @author: Jeffrey Chan, 2013
%

    assert(length(vPmf1) == length(vPmf2));
    assert(sum(vPmf1) == 1);
    assert(sum(vPmf2) == 1);
    
    vNz1 = find(vPmf1);
    vNz2 = find(vPmf2);

    div = sum(vPmf1(vNz1) .* log2(vPmf1(vNz1) ./ (0.5.*(vPmf1(vNz1) + vPmf2(vNz1))))) + sum(vPmf2(vNz2) .* log2(vPmf2(vNz2) ./ (0.5.*(vPmf1(vNz2) + vPmf2(vNz2)))));

end % end of fuction

