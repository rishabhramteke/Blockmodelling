function plotSynResult(sBaseDir, sMatWildcard, sXlabel, sYLabel, xRange)
%
% Plots one plot per data.
%
% mData - 3xn matrix, where first column is x-values, 2nd is mean and 3rd
% is the variance.
% sXLabel - x axis label
% sYLabel - y axis label
%


% get the list of snapshot/matrices filenames
stMatFilenames = dir(fullfile(sBaseDir, sMatWildcard));

stMatFilenames.name
% load them one by ''one
for i = 1 : size(stMatFilenames,1)
    mData = csvread(fullfile(sBaseDir, stMatFilenames(i).name), 0, 0);
    
    figure;
    errorbar(mData(:,1), mData(:,2), mData(:,3), '-b');
    set(gca, 'FontSize', 40);
    set(gca, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 2);
    set(gca, 'XLim', xRange);
    ylabel(sYLabel);
    xlabel(sXlabel);    
    % save file
    baseIndex = regexp(stMatFilenames(i).name, '\.');
    filename = stMatFilenames(i).name(1:baseIndex-1);
    saveas(gca, strcat(filename, '.fig'), 'fig');
    saveas(gca, strcat(filename, '.jpg'), 'jpg');
end








end