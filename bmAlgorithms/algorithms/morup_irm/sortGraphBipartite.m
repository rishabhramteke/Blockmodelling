function [A,Z1,Z2,eta,perm1,perm2]=sortGraphBipartite(A,Z1,Z2,eta)
% Function to sort graph accourding to clustering structure derived from
% bi-partite IRM model. The clusters are sorted according to the cluster sizes.
%
% Usage:
%   [A,Z1,Z2,eta,perm1,perm2]=sortGraphBipartite(A,Z1,Z2,eta)
%
% Input:
%    A      Adjacency matrix or multigraph
%    Z1     Clustering assignment matrix of mode 1
%    Z2     Clustering assignment matrix of mode 2
%    eta    matrix of between cluster relations
%
% Output: 
%    A      sorted Adjacency matrix or multigraph
%    Z1     sorted clustering assignment matrix of mode 1
%    Z2     sorted clustering assignment matrix of mode 2
%    eta    sorted matrix of between cluster relations
%    perm1  the permutation of the indices of mode 1
%    perm2  the permutation of the indices of mode 1
%
% Written by Morten Mørup


% Sort Z1, Z2 and eta according to cluster size
sumZ1=sum(Z1,2);
[val,ind1]=sort(sumZ1,'descend');
Z1=Z1(ind1,:);
sumZ2=sum(Z2,2);
[val,ind2]=sort(sumZ2,'descend');
Z2=Z2(ind2,:);
eta=eta(ind1,:);
eta=eta(:,ind2);

% Sort A and Z1 and Z2 according to cluster assignment
noc1=size(Z1,1);
[val,perm1]=sort((2.^(0:noc1-1))*Z1,'ascend');
Z1=Z1(:,perm1);
if iscell(A)
    for n=1:length(A)
        A{n}=A{n}(perm1,:);
    end
else
     A=A(perm1,:);
end
noc2=size(Z2,1);
[val,perm2]=sort((2.^(0:noc2-1))*Z2,'ascend');
Z2=Z2(:,perm2);
if iscell(A)
    for n=1:length(A)
        A{n}=A{n}(:,perm2);
    end
else
     A=A(:,perm2);
end



