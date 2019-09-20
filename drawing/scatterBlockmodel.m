function scatterBlockmodel(vX, vY, vVertLineLocs, vHorLineLocs)
%
% Draws the blockmodel adjacency matrixvVertLineLocs
% mAdj - adjacency matrix, ordered into blocks already
% vVertLineLocs - location of vertical lines
% vHorLineLocs - location of horizonal lines
%


scatter(vX, vY, 1, 'filled');

xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');


% verical line
if size(vVertLineLocs,2) > 0
    line([vVertLineLocs; vVertLineLocs], xlim, 'LineStyle', '--', 'Color', 'red');
end
if size(vHorLineLocs,2) > 0
    line(ylim, [vHorLineLocs; vHorLineLocs], 'LineStyle', '--', 'Color', 'red');
end

end % end of function