function [dist] = compareOverlapBM(sPosFile1, sGraphFile1, sPosFile2, sGraphFile2,...
    bHeader1, bHeader2, bListFormat, vNum1, vNum2, bMatlabFormat, bSparse, bAddOneToIndex, bAddOneToPosVerts, sReconDist, sPosDist, sCombDist, sOutFile, combParam)
%
% Computes the positional and density distances between clustering/bm 1 and
% 2.
%
% sReconDist - Reconstruction distance to use (if not used, put '').
% sPosDist - Position distance to use (if not used, put '').
% sCombDist - Combining distance to use (if not used, put '').
%

% load positions
[mVertPosMap1, aPartInfo1] = loadPositions(sPosFile1, bHeader1, bListFormat, vNum1, bAddOneToPosVerts);
[mVertPosMap2, aPartInfo2] = loadPositions(sPosFile2, bHeader2, bListFormat, vNum2, bAddOneToPosVerts);


% load graph file
[mSparseGraph1] = loadMatlabGraph(sGraphFile1, vNum1, vNum1, bMatlabFormat, bSparse, bAddOneToIndex);
[mSparseGraph2] = loadMatlabGraph(sGraphFile2, vNum2, vNum2, bMatlabFormat, bSparse, bAddOneToIndex);


dist = compareOverlapBMNoFile(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sReconDist, sPosDist, sCombDist, sOutFile, combParam);


% if (size(sReconDist,2) == 0 && size(sPosDist,2) == 0)
%     error('No distance measure specified');
% else
% 
%     % open approprite file for writing
%     fOut = fopen(sOutFile,'w');
%     
%     % combine distance option
%     if (size(sCombDist,2) > 0)
%         assert(size(sReconDist,2) > 0 && size(sPosDist,2) > 0);
%         denDis = bmOverlapDenDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sReconDist);
%         posDis = bmOverlapPosDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sPosDist);
%         combDis = bmCompareCombDist(posDis, denDis, sCombDist, combParam);
%         fprintf(fOut,'%.12f\n', combDis);
%         dist = combDis;
%     else
%         if (size(sReconDist,2) > 0)
%             denDis = bmOverlapDenDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sReconDist);
%         end
% 
%         if (size(sPosDist,2) > 0)
%             posDis = bmOverlapPosDist(mVertPosMap1, mSparseGraph1, mVertPosMap2, mSparseGraph2, sPosDist);
%         end
%         
%         % write out to file
%         if (size(sReconDist,2) > 0 && size(sPosDist,2) > 0)
%             fprintf(fOut,'%.12f %.12f\n', posDis, denDis);
%             dist = [posDis, denDis];
%         elseif  (size(sReconDist,2) > 0)
%             fprintf(fOut,'%.12f \n', denDis);
%             dist = denDis;
%         elseif (size(sPosDist,2) > 0)
%             fprintf(fOut,'%.12f \n', posDis);
%             dist = posDis;
%         end        
%     end
% 
%     fclose(fOut);
% end

% need to exit if run as script
exit;

end % end of function




