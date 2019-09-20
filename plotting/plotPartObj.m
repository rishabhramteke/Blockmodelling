function plotPartObj(bTrueResult, sTitle, varargin)
%
% Plot the input data as subplots.
%
% bTrueResult:  boolean to indicate if the dataset loaded has true
% distribution/partition results.
%

optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin < 0
    disp(' plotPartObj(varargin)');
    exit;
end


% plot
figure;

% calculate the size of the subplot grid (2 subplots per row)
totalRow = ceil(optargin / 2);
if bTrueResult
    totalRow = ceil(optargin);
end
%totalCol = 1;

    totalCol = 2;


currRow = 0;
% loop through the data
for i = 1 : optargin       
    mData = varargin{i};
    mMin = mData(:,2:6);
    mMax = mData(:,7:11);
    mMean = mData(:,12:16);
    mStd = mData(:,17:21);
    mMedian = mData(:,22:26);
    if bTrueResult
        mMin = mData(:,2:7);
        mMax = mData(:,8:13);
        mMean = mData(:,14:19);
        mStd = mData(:,20:25);
        mMedian = mData(:,26:31);
        trueTotalCost = mData(:,32);
    end

    
    
    
    subplot(totalRow, totalCol, currRow + i);
    
        
    if bTrueResult
        plot(mData(:,1), mMin(:,1), '-r',  mData(:,1), mMin(:,2), '-b', mData(:,1), mMin(:,3) - mMin(:,4), '--y', mData(:,1), mMin(:,4), '-.g', mData(:,1), trueTotalCost, ':b', mData(:,1), mMean(:,6), ':b');
        legend('total cost', 'dist cost', 'block cost', 'block entropy', 'true cost', 'avg(random partitions)');
    else
        plot(mData(:,1), mMin(:,1), '-r',  mData(:,1), mMin(:,2), '-b', mData(:,1), mMin(:,3) - mMin(:,4), '--y', mData(:,1), mMin(:,4), '-.g', mData(:,1), mMean(:,5), ':b');
        legend('total cost', 'dist cost', 'block cost', 'block entropy', 'avg(random partitions)');
    end
    
    ylabel('Objective cost (bits)');
    xlabel('Number of tested partitions');
    % 'Best objective costs, True partition = 10'
    title(sTitle);    
    
    
    if bTrueResult
        subplot(totalRow, totalCol, currRow + i+1);
        [AX,H1,H2] = plotyy(mData(:,1), mMin(:,1), mData(:,1), mMin(:,5));
        set(get(AX(1),'Ylabel'),'String','Objective cost (bits)');
        set(get(AX(2),'Ylabel'),'String','Variation of Information (bits)');
        set(H1,'LineStyle','-');
        set(H2,'LineStyle','--');    
        xlabel('Number of tested partitions');
        title('Best objective cost + clustering distance, True partition = 10');
        legend('total cost', 'variation of information');
    end

    if bTrueResult
        currRow = currRow + 1;
    end
end


end % end of function