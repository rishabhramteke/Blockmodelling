function plotDendrogramsSummerSchool(cmDistances, sLinkageType)
%
% Plots one dendrogram per cell in cmDistances (i.e., one dendrogram per
% each similarity matrix done for a measure.
%

for i = 1 : size(cmDistances, 2)
    Z = linkage(squareform(cmDistances{i} + cmDistances{i}', 'tovector'), sLinkageType);
    figure;
    sTitle = sprintf('%d', i);
    title(sTitle);
    H = dendrogram(Z, 'labels', {'Role', 'SQB', 'IRM', 'SBM3', 'SBM4', 'SBM8'});
    set(H,'LineWidth',3)
end