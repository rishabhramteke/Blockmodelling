function div =  klDivergence(vPmf1, vPmf2, vWeightDistValues)
%
% Computes the KL divergence between the two pmfs.
%
% @author: Jeffrey Chan, 2013
%

    assert(length(vPmf1) == length(vPmf2));
    assert(abs(sum(vPmf1)- 1) <= 10 * eps);  
    assert(abs(sum(vPmf2) -1) <= 10 * eps);
    
    vPmf1Alt = vPmf1 + eps;
    vPmf2Alt = vPmf2 + eps;
    
    div = sum(vPmf1Alt .* log2(vPmf1Alt ./ vPmf2Alt));

end % end of fuction