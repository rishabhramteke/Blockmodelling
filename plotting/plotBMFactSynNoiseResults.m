function plotBMFactSynNoiseResults(sBaseDir, xCol, yCol, yError, sXlabel, sYLabel, bLogY, vXRange)
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

        
    %hard + comparison with existing
%     ccAlgorNames = {...
%         {'coordDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'projGradDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'bolongM','hardCIncr','euclidean'},...
%         {'saHamilM', 'saHamilM', 'reichEqual'}
%         };
% %    
% %             
%     vsLegendNames = {...
%         'Coor-H-Ad',...
%         'Grad-H-Ad',...
%         'Long-H',...
%         'Reic-H'
%         };
    
    % objective
%     ccAlgorNames = {...
%         {'projGradDescM','hardCIncr','euclidean'},...
%         {'coordDescExactM','hardCIncr','euclidean'},...
%         {'projGradDescM','hardCIncr','bmEuclidean'},...
%         {'coordDescM','hardCIncr','bmEuclidean'},...
%         {'projGradDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'coordDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...        
%         {'projGradDescM','hardCIncrAdjEqual','bmEuclideanAdjEqual'},...
%         {'coordDescM','hardCIncrAdjEqual','bmEuclideanAdjEqual'}
%         };
%

%     ccAlgorNames = {...
%         {'projGradDescM','hardCIncrAdjEqual','bmEuclideanAdjEqual'},...
%         {'projGradDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'projGradDescM','hardCIncr','bmEuclidean'},...
%         {'projGradDescM','hardCIncr','euclidean'}
%         };
%

%     ccAlgorNames = {...
%         {'coordDescBMM','hardCIncrAdjEqual','bmEuclideanAdjEqual'},...
%         {'coordDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'coordDescBMM','hardCIncr','bmEuclidean'},...
%         {'coordDescExactM','hardCIncr','euclidean'}
%         };

% 
%     vsLegendNames = {...
%         'Grad-H',...
%         'Coor-H',...
%         'Grad-H-Con',...
%         'Coor-H-Con',...
%         'Grad-H-Adj',...
%         'Coor-H-Adj',...
%         'Grad-H-AdjCon',...
%         'Coor-H-AdjCon'
%         };

%     vsLegendNames = {...
%         'Grad-H-AdCn',...
%         'Grad-H-Ad',...
%         'Grad-H-Cn',...
%         'Grad-H',...
%         };
% %
%     vsLegendNames = {...
%         'Coor-H-AdCn',...
%         'Coor-H-Ad',...
%         'Coor-H-Cn',...
%         'Coor-H'
%         };


% 

%     % soft
%     ccAlgorNames = {...
%         {'bolongM','softCMultLong','euclidean'},...
%         {'dingM','softCMultDing','euclidean'},...
%         {'coordDescM','softCCoordDesc','euclidean'},...
%         {'projGradDescM','softCGradDesc','euclidean'},...
%         {'coordDescM','softCCoordDescAdjEqual','euclideanAdjEqual'},...
%         {'projGradDescM','softCGradDescAdjEqual','euclideanAdjEqual'},...
%         {'coordDescBMM','softCGradDescAdjEqual','bmEuclideanAdjEqual'},...
%         {'projGradDescM','softCGradDescAdjEqual','bmEuclideanAdjEqual'}
%         };
% 
    ccAlgorNames = {...
        {'coordDescM','softCCoordDescAdjEqual','euclideanAdjEqual'},...
        {'projGradDescM','softCGradDescAdjEqual','euclideanAdjEqual'},...
        {'bnmtfM','softCBnmtf','euclidean'},...
        {'bolongM','softCMultLong','euclidean'},...
        {'dingM','softCMultDing','euclidean'},...
        };
%             
% 
%     vsLegendNames = {...
%         'BNMTF',...
%         'Long-S',...
%         'Ding-S',...
%         'Coor-S',...
%         'Grad-S',...
%         'Coor-S-Adj',...
%         'Grad-S-Adj',...
%         'Coor-S-AdjCon',...
%         'Grad-S-AdjCon'
%         };


    vsLegendNames = {...
        'Coor-S-Ad',...
        'Grad-S-Ad',...
        'BNMTF',...
        'Long-S',...
        'Ding-S'
        };



    % graph size
%     ccAlgorNames = {...
%         {'coordDescExactM','hardCIncr','euclidean'},...
%         {'projGradDescM','hardCIncr','euclidean'},...
%         {'coordDescM','hardCIncr','bmEuclidean'},...
%         {'projGradDescM','hardCIncr','bmEuclidean'},...
%         {'coordDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'projGradDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'coordDescM','hardCIncrAdjEqual','bmEuclideanAdjEqual'},...
%         {'projGradDescM','hardCIncrAdjEqual','bmEuclideanAdjEqual'},...
%         };
%     
%     vsLegendNames = {...
%         'Coor-H',...
%         'Grad-H',...
%         'Coor-H-Cn',...
%         'Grad-H-Cn',...
%         'Coor-H-Ad',...
%         'Grad-H-Ad',...
%         'Coor-H-AdCn',...
%         'Grad-H-AdCn',...
%         };    
    
    
    
%     ccAlgorNames = {...
%         {'bolongM','hardCIncr','euclidean'},...
%         {'coordDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'projGradDescM','hardCIncrAdjEqual','euclideanAdjEqual'},...
%         {'bnmtfM','softCBnmtf','euclidean'},...
%         {'bolongM','softCMultLong','euclidean'},...
%         {'dingM','softCMultDing','euclidean'},...
%         };
% 
%             
% 
%     vsLegendNames = {...
%         'Long-H',...
%         'Coor-H-Ad',...
%         'Grad-H-Ad',...
%         'BNMTF',...
%         'Long-S',...
%         'Ding-S',...
%         };




    figure;
    hold on;

    cc=hsv(length(ccAlgorNames));
    
    vLineMarker = {'+', 'o', '*', '^', 'x', 's', 'd', '.', 'v', '>', '<', 'p', 'h'};
%     vLineStyle = {'-', '--', ':', '-.'};
    vLineStyle = {'-', '--', '-.'};
        
    % default line settings
    lineWidth = 3;
    markerSize = 3;
    fontSize = 24;


%     vsLegendNames = cell(1, length(ccAlgorNames));
    % load them one by ''one
    for i = 1 : length(ccAlgorNames)
        cNames = ccAlgorNames{i}
        % get the filename
        stMatFilenames = dir(fullfile(sBaseDir, strcat('*', cNames{1}, '_', cNames{2}, '_', cNames{3}, '*')))
        assert(length(stMatFilenames) == 1);
        
%         vStart = regexp(stMatFilenames(1).name, '_');
%         sName = strcat(stMatFilenames(1).name(vStart(1)+1:vStart(2)-1), '-',...  
%             stMatFilenames(1).name(vStart(2)+1:vStart(3)-1), '-',...
%             stMatFilenames(1).name(vStart(3)+1:vStart(4)-1));
%         vsLegendNames{i} = sName;
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
        
        errorbar(mData(:,xCol), mData(:,yCol), mData(:,yError), 'color',cc(i,:), 'marker', vLineMarker{m}, 'LineStyle', vLineStyle{s}, 'LineWidth', lineWidth, 'markers', markerSize);

%         set(gca, 'XLim', xRange);
  

%         saveas(gca, strcat(filename, '.fig'), 'fig');
%         saveas(gca, strcat(filename, '.jpg'), 'jpg');
    end

    set(gca, 'FontSize', fontSize);
    set(gca, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 4);
    if bLogY
        set(gca, 'YScale', 'log');
    end
    ylabel(gca, sYLabel);
    xlabel(gca, sXlabel);  
    xlim(vXRange);
    legend(vsLegendNames);
    hold off;
    
end