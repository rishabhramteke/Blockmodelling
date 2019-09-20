function quadLossNMFTesting(mMlData, mMlItem, vKRange)
%
% @author: Jeffrey Chan, 2014
%

    for k = vKRange(1) : vKRange(2)
        [mW, mH] = nnmf(mMlData, k);
   
        plotInfo(mMlData, mMlItem, mW, mH);
    end
    

end % end of function




function plotInfo(mMlData, mMlItem, mW, mH)
%
% Plots various stats about each algorithm results over the iterations of their
% runs.  One plot per run.
%

        
    % colour of lines
    cColour = {'b', 'r', 'y', 'g', 'm', 'c', 'k', 'b', 'r', 'y', 'g', 'm', 'c', 'k', 'b', 'r', 'y', 'g', 'm', 'c', 'k'};
    % plot the various things
    rowNum = 2;
    colNum = 3;
    
    mDiff = (mMlData - mW * mH).^2;
    
    figure;
    figIndex = 1;
    
    % heatmap
    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    imagesc(mDiff);
    title('Headmap: Squared Error');
    colormap(flipud(gray));
    
    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    hist(mDiff(:), 100);   
    title('Histogram: Squared Error');
    set(gca, 'YScale', 'log');
    
    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    imagesc(mMlData);  
    title('Original Data');
    colormap(flipud(gray));
    
    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    imagesc(mW * mH);  
    title('W * H');
    colormap(flipud(gray));
    
    

    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    imagesc(mMlItem(:,2:end));  
    title('Item categories');
    colormap(flipud(gray));    
    
    
    
    subplot(rowNum, colNum, figIndex);
    figIndex = figIndex + 1;
    vData = zeros(size(mDiff,1) * size(mDiff,2) * length(mMlItem(1,2:end)), 1);
    vTag = zeros(size(mDiff,1) * size(mDiff,2) * length(mMlItem(1,2:end)), 1);
    % square error per tag
    currIndex = 1;
    for k = 2 : size(mMlItem(:,2:end),2)
        vI = mMlItem(:,k) > 0;
        mTemp = mDiff(vI, :);
        vTemp = mTemp(:);
        vData(currIndex : currIndex + length(vTemp)-1) = vTemp;
        vTag(currIndex : currIndex + length(vTemp)-1) = k;
        currIndex = currIndex + length(vTemp);
    end
    % get rid of tail
    vData(currIndex+1:end) = [];
    vTag(currIndex+1:end) = [];
    boxplot(vData, vTag, 'plotstyle', 'traditional');
    
    
    
    
    
end % end of function
