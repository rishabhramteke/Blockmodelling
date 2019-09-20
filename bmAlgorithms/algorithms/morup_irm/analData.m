type='UnDirected';  % Graph generated will be undirected
J=500;              % Size of graph
alpha=5;            % parameter for the CRP
bp=[1 1];           % prior pseudo-link counts within and between communities
bn=[1 1];           % prior pseudo-nonlink counts within and between communities

% Generate graph according to the IRM model with the specified parameters
[A,Z_true,eta_true,A_sorted,Z_sorted,eta_sorted,perm]=generateGraphCRPUnipartite(J,alpha,bp,bn,type);

% Define entries in graphs to be treated as missing in the inference
pct_missing=2.5;
[W,class]=createValidationData(A,pct_missing,type);

% Infer the parameters of the IRM model from the generated graph
noc=5; % Initial number of components
opts.type=type;
[L,cpu_time,Z_estimated,eta_estimated,sample,West]=IRMUnipartite(A,W,noc,opts);

% Plot the results
plotSyntheticResults(A,West,Z_true,eta_true,Z_estimated,eta_estimated);