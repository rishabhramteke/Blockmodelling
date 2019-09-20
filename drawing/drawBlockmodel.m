function drawBlockmodel(mAdj, vVertLineLocs, vHorLineLocs, varargin)
%
% Draws the blockmodel adjacency matrixvVertLineLocs
% mAdj - adjacency matrix, ordered into blocks already
% vVertLineLocs - location of vertical lines
% vHorLineLocs - location of horizonal lines
%
% @author: Jeffrey Chan, 2012
%

    mNewAdj = mat2gray(mAdj);

    imagesc(mNewAdj);
    colormap(flipud(gray));
    caxis([0,1]);

    xlim = get(gca, 'xlim');
    ylim = get(gca, 'ylim');

    if length(varargin) == 2
        set(gca,'XTickLabel', varargin{1});
        set(gca,'YTickLabel', varargin{2});
        set(gca,'XTick', 1:length(varargin{1}));
        set(gca,'YTick', 1:length(varargin{2}));
    end


    % verical line
    if size(vVertLineLocs,2) > 0
        line([vVertLineLocs; vVertLineLocs], xlim, 'LineStyle', '--', 'Color', 'red', 'LineWidth', 3);
    end
    if size(vHorLineLocs,2) > 0
        line(ylim, [vHorLineLocs; vHorLineLocs], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 3),;
    end

end % end of function