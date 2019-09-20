function plotBMFactSynSparsityResults(sBaseDir, xCol, yCol, yError, sXlabel, sYLabel, bLogY)
%
% Plots one plot per data.
%
% mData - 3xn matrix, where first column is x-values, 2nd is mean and 3rd
% is the variance.
% sXLabel - x axis label
% sYLabel - y axis label
%


    % get the list of snapshot/matrices filenames
%     stMatFilenames = dir(fullfile(sBaseDir, sResultsWildcard));


    ccAlgorNames = {...
        {'bnmtfM','softCBnmtf','euclidean'},...
        {'bolongM','softCMultLong','euclidean'},...
        {'dingM','softCMultDing','euclidean'},...
        {'coordDescM','softCCoordDescAdjEqual','bmEuclideanAdjEqual'},...
        {'projGradDescM','softCGradDescAdjEqual','bmEuclideanAdjEqual'}
        };

    figure;
    hold on;

    cc=hsv(length(ccAlgorNames));
    
    vLineMarker = {'+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
    vLineStyle = {'-', '--', ':', '-.'};
        



    vsLegendNames = cell(1, length(ccAlgorNames));
    % load them one by ''one
    for i = 1 : length(ccAlgorNames)
        cNames = ccAlgorNames{i};
        % get the filename
        stMatFilenames = dir(fullfile(sBaseDir, strcat('*', cNames{1}, '_', cNames{2}, '_', cNames{3}, '*')));
        assert(length(stMatFilenames) == 1);
        
        vStart = regexp(stMatFilenames(1).name, '_');
        sName = strcat(stMatFilenames(1).name(vStart(1)+1:vStart(2)-1), '-',...  
            stMatFilenames(1).name(vStart(2)+1:vStart(3)-1), '-',...
            stMatFilenames(1).name(vStart(3)+1:vStart(4)-1));
        vsLegendNames{i} = sName;
        mData = csvread(fullfile(sBaseDir, stMatFilenames(1).name), 0, 0);
    
        % determine marker to use
        m = mod(i, length(vLineMarker));
        if i / length(vLineMarker) >= 1
            m = m + 1;
        end
        
        s = mod(i, length(vLineStyle));
        if i / length(vLineStyle) >= 1
            s = s + 1;
        end 
        
        errorbar(mData(:,xCol), mData(:,yCol), mData(:,yError), 'color',cc(i,:), 'marker', vLineMarker{m}, 'LineStyle', vLineStyle{s});

%         set(gca, 'XLim', xRange);
  

%         saveas(gca, strcat(filename, '.fig'), 'fig');
%         saveas(gca, strcat(filename, '.jpg'), 'jpg');
    end

    set(gca, 'FontSize', 40);
    set(gca, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 2);
    if bLogY
        set(gca, 'YScale', 'log');
    end
    ylabel(gca, sYLabel);
    xlabel(gca, sXlabel);  
    xlim([-0.05 1.0]);
    legend(vsLegendNames);
    hold off;
    
end