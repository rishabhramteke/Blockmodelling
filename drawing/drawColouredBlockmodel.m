function drawColouredBlockmodel(mAdj, mMembership, mPosColour, vVertLineLocs, vHorLineLocs, varargin)
%
% Draws the blockmodel as an block pixel image, with different colours for each
% position and block.  Each position is assigned the colour specified in
% mPosColour and then the colour of the block is the resultant mixture of the
% two.
%
% mAdj - adjacency matrix, ordered into blocks already
% mPosColour - the RGB representation of the colour to use for each position.
% vVertLineLocs - location of vertical lines
% vHorLineLocs - location of horizonal lines
% varargin - xlabels and ylabels.
%


    % construct colour 3d matrix 
    m3ColourAdj = ones(size(mAdj,1), size(mAdj,2), 3);
        
    [vNonZeroRow, vNonZeroCol] = find(mAdj > 0);
    
    % colour each edge in mAdj
    for i = 1 : length(vNonZeroRow)
        % find position
        rowPos = mMembership(vNonZeroRow(i),:) > 0;
        colPos = mMembership(vNonZeroCol(i),:) > 0;
        vCombinedColour = (mPosColour(rowPos,:) + mPosColour(colPos,:)) / 2;
        m3ColourAdj(vNonZeroRow(i), vNonZeroCol(i), :) = cat(3, vCombinedColour(1), vCombinedColour(2), vCombinedColour(3));
    end

    
    image(m3ColourAdj);


    xlim = get(gca, 'xlim');
    ylim = get(gca, 'ylim');
    
    if length(varargin) == 2
        set(gca,'XTickLabel', varargin{1});
        set(gca,'YTickLabel', varargin{2});
        set(gca,'XTick', 1:length(varargin{1}));
        set(gca,'YTick', 1:length(varargin{2}));
    end

    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');

    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

    
    
    % vertical line
    if size(vVertLineLocs,2) > 0
        line([vVertLineLocs; vVertLineLocs], xlim, 'LineStyle', '--', 'Color', 'red', 'LineWidth', 3);
    end
    if size(vHorLineLocs,2) > 0
        line(ylim, [vHorLineLocs; vHorLineLocs], 'LineStyle', '--', 'Color', 'red', 'LineWidth', 3);
    end

end % end of function