function plotSyntheticResults(A,West,Z_true,eta_true,Z_estimated,eta_estimated)
% function to visualize the analysis of synthetically generated IRM data
%
% Usage:
%    plotSyntheticResults(A,Ztrue,Zestimated)
%
% Input:
%   A            Generated adjancency matrix
%   West         Link-predicted values
%   Ztrue        True assignment matrix
%   Zestimated   Estimated assignment matrix
%
% Written by Morten Mørup

J=size(A,1);
figure;
subplot(2,3,1);
mySpyPlot(A,1000/J);
title('Generated Graph','FontWeight','Bold')
subplot(2,3,2);
[val,ind]=sort(sum(Z_true,2),'descend');
Z_true=Z_true(ind,:);
eta_true=eta_true(ind,:);
eta_true=eta_true(:,ind);
[A_sorted,Z_sorted,eta_sorted]=sortGraphUnipartite(A,Z_true,eta_true);
mySpyPlot(A_sorted,1000/J,Z_sorted,Z_sorted,eta_sorted);
title('Correctly Sorted Generated Graph','FontWeight','Bold')
subplot(2,3,3);
[A_sorted,Z_sorted,eta_sorted]=sortGraphUnipartite(A,Z_estimated,eta_estimated);
mySpyPlot(A_sorted,1000/J,Z_sorted,Z_sorted,eta_sorted);
title('Estimated Sorting of the Generated Graph','FontWeight','Bold')
subplot(2,3,4);
imagesc(-Z_true*Z_estimated');  axis equal; axis tight; title(['Ztrue*Zestimated^T, NMI=' num2str(round(calcNMI(Z_true,Z_estimated)*100)/100) ', AUC=' num2str(round(calcAUC(West,A)*100)/100)],'FontWeight','Bold')
subplot(2,3,5);
imagesc(-Z_true);  axis off; title('True Assigment Matrix Z','FontWeight','Bold')
subplot(2,3,6);
imagesc(-Z_estimated); axis off; title('Estimated Assignment Matrix Z','FontWeight','Bold')
