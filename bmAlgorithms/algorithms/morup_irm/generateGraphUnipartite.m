function [A,Z,eta,mu,A_sorted,Z_sorted,eta_sorted,mu_sorted,perm]=generateGraphUnipartite(noc,J,alpha,bp,bn,type)
% Script to generate unipartite graphs according to the Finite Relational
% Model (i.e. Stochastic Block Model)
%
% Usage:
%   [A,Z,eta,mu,A_sorted,Z_sorted,eta_sorted,mu_sorted,perm]=GenerateGraphUnipartite(noc,J,alpha,bp,bn,type)
% 
% Input:
%   noc         number of components
%   J           size of graph 
%   alpha       parameter for the Chinese Restaurant Process (CRP)
%   bp          1x2 vector (default [10 1]) where 
%               bp(1): indicate within community prior link count
%               bp(2): indicate between community prior link count
%   bn          1x2 vector (default [1 10]) where 
%               bn(1): indicate within community prior non-link count
%               bn(2): indicate between community prior non-link count
%   type        'UnDirected' (default) or 'Directed'
%
% Output:
%   A           Generated graph that has not been sorted, i.e. A=Bernoulli(Z'*eta*Z)
%   Z           noc x J generated assignment matrix, i.e. Z~Discrete(mu)
%   eta         noc x noc generated group relations, i.e. eta_{lm}~Beta(Bp,Bn)
%   mu          noc x 1 cluster probability vector, i.e. mu~Dirichlet(alpha)
%   A_sorted    sorted J x J Adjacency matrix of the generated graph
%   Z_sorted    sorted noc x J generated assignment matrix
%   eta_sorted  sorted noc x noc generated group relations
%   mu          sorted noc x 1 cluster probability vector
%   perm        permutation matrix, i.e. A_sorted=A(perm,perm)
%
% Written by Morten Mørup

if nargin<6
    type='UnDirected';
end
if nargin<5
    bn=[1 10];    
end
if nargin<4
    bp=[10 1];
end
if nargin<3
    alpha=1;
end

if length(alpha)==1
    alpha=alpha*ones(noc,1); % Parameter for dirichlet distribution imposed on mu
end

%%%%%%%%%%%%%%%%%%
% Generate graph %
%%%%%%%%%%%%%%%%%%
% mu~dirichlet(alpha)
mu=dirichletrnd(alpha);

% z_i ~ Discrete(mu)
pp=cumsum(mu);
Z=zeros(noc,J);
for i=1:J
   ind=find(rand<pp,1,'first');
   Z(ind,i)=1;
end

% eta_lm ~ Beta(bp,bn)
eta=betarnd(bp(1)*eye(noc)+bp(2)*(ones(noc)-eye(noc)),bn(1)*eye(noc)+bn(2)*(ones(noc)-eye(noc)));
if strcmp(type,'UnDirected') % Force eta=eta' for UnDirected graphs
    eta=triu(eta,1)+triu(eta,1)'+diag(diag(eta));     
end

% A_ij ~ Bernoulli(z_i'*eta*z_j)
A=sparse((Z'*eta*Z)>rand(J));
if strcmp(type,'UnDirected')
    A=triu(A,1);
    A=A+A';
else
    A=A-diag(diag(A));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display generated graph %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1);
mySpyPlot(A,1000/J);
title('Generated Graph','FontWeight','Bold')
subplot(2,2,2);
[A_sorted,Z_sorted,eta_sorted,perm,mu_sorted]=sortGraphUnipartite(A,Z,eta,mu);
mySpyPlot(A_sorted,1000/J,Z_sorted,Z_sorted,eta_sorted);
title('Sorted Generated Graph','FontWeight','Bold')
subplot(2,2,3);
imagesc(-Z_sorted);  axis off; title('Generated Sorted Clustering Assigment Matrix Z','FontWeight','Bold')
subplot(2,2,4);
imagesc(-eta_sorted); axis equal; axis tight; title('Generated Between Cluster Relations \eta','FontWeight','Bold')


