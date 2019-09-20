function [dist] = compareOverlapBMNoFile(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2,...
    sReconDist, sPosDist, sCombDist, sOutFile, combParam)
%
% Computes the positional and density distances between clustering/bm 1 and
% 2.
%
% sReconDist - Reconstruction distance to use (if not used, put '').
% sPosDist - Position distance to use (if not used, put '').
% sCombDist - Combining distance to use (if not used, put '').
%



if (size(sReconDist,2) == 0 && size(sPosDist,2) == 0)
    error('No distance measure specified');
else

    % open approprite file for writing
    fOut = fopen(sOutFile,'w');
    
    % combine distance option
    if (size(sCombDist,2) > 0)
        assert(size(sReconDist,2) > 0 && size(sPosDist,2) > 0);
        denDis = bmOverlapDenDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sReconDist);
        posDis = bmOverlapPosDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sPosDist);
        combDis = bmCompareCombDist(posDis, denDis, sCombDist, combParam);
        fprintf(fOut,'%.12f\n', combDis);
        dist = combDis;
    else
        if (size(sReconDist,2) > 0)
            denDis = bmOverlapDenDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sReconDist);
        end

        if (size(sPosDist,2) > 0)
            posDis = bmOverlapPosDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sPosDist);
        end
        
        % write out to file
        if (size(sReconDist,2) > 0 && size(sPosDist,2) > 0)
            fprintf(fOut,'%.12f %.12f\n', posDis, denDis);
            dist = [posDis, denDis];
        elseif  (size(sReconDist,2) > 0)
            fprintf(fOut,'%.12f \n', denDis);
            dist = denDis;
        elseif (size(sPosDist,2) > 0)
            fprintf(fOut,'%.12f \n', posDis);
            dist = posDis;
        end        
    end

    fclose(fOut);
end


end % end of function




