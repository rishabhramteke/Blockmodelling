function plotComparisonResults(cmResults, vIndices, vNames, sXlabel, sYLabel, bYLog, vXRange)
%
% INPUT:
% cmResults     - cell of matrices of results.
% vIndices      - vector of indices to use.
% vNames        - vector of names.
%
    
    figure;
    hold on;
    
    cc=hsv(length(vIndices));
    
    vLineMarker = {'+', 'o', '*', '^', 'x', 's', 'd', '.', 'v', '>', '<', 'p', 'h'};
%     vLineStyle = {'-', '--', ':', '-.'};
    vLineStyle = {'-', '--', '-.'};
        
    % default line settings
    lineWidth = 3;
    markerSize = 3;
    fontSize = 24;    
    
    
    % number of positions moved
    xSize = size(cmResults{1},2);
    
    for i = 1 : length(vIndices)
        % determine marker to use
        m = mod(i, length(vLineMarker));
        if i / length(vLineMarker) >= 1
            m = m + 1;
        end
        
        % determine line style to use
        s = mod(i, length(vLineStyle));
        if i / length(vLineStyle) >= 1
            s = s + 1;
        end 
        
        errorbar(1:xSize, mean(cmResults{i},1), std(cmResults{i},0,1), 'color',cc(i,:), 'marker', vLineMarker{m}, 'LineStyle', vLineStyle{s}, 'LineWidth', lineWidth, 'markers', markerSize);
    end
    
    set(gca, 'FontSize', fontSize);
    set(gca, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 4);
    if bYLog
        set(gca, 'YScale', 'log');
    end
    ylabel(gca, sYLabel);
    xlabel(gca, sXlabel);  
    xlim(vXRange);
    legend(vNames);    
    hold off;

end % end of function





   