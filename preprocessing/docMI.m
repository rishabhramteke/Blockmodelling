function [] = docMI(mDocTermFreq, vDocClass, classNum)
%
% Computes the pairwise mutual information.
%
% mDocTermFreq -    The matrix of document-term-frequency.
% vDocClass -       The vector of document-class label.
% classNum -        Class number.
%





% perform a join between mDocTermFreq and vDocClass to produce class-term-freq

% get all the terms together
mDocTermFreqSort = sortrows(mDocTermFreq, 2);


for d = 1 : length(vDocClass)
    
end

    

end % end of function