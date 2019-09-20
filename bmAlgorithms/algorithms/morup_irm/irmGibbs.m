function Z = irmGibbs(X,T,a,b,A)

% Simple script for sampling z according to the IRM model based on
% UnDirected graphs

N = size(X,1); z = true(N,1); Z = cell(T,1);                                % Initialization

for t = 1:T                                                                 % For each Gibbs sweep        
    for n = 1:N                                                             % For each node in the graph
        z(n,:) = false;                                                     % Disregard node n 
        K = size(z,2);                                                      % No. of components
        m = sum(z)'; M = repmat(m,1,K);                                     % No. of nodes in each component
        M1 = z'*X*z-...                                                     % No. of links between components        
            diag(sum(X*z.*z)/2);                         
        M0 = m*m'-diag(m.*(m+1)/2)-M1;                                      % No. of non-links between components
        r = z'*X(:,n); R = repmat(r,1,K);                                   % No. of links from node n
        logP = sum([betaln(M1+R+a,M0+M-R+b)-betaln(M1+a,M0+b) ...
            betaln(r+a,m-r+b)-betaln(a,b)],1)' + log([m; A]);               % Log probability of n belonging to  existing and new component                           
        P = exp(logP-max(logP));                                            % Convert from log probability
        i = find(rand<cumsum(P)/sum(P),1);                                  % Random component according to P
        z(n,i) = true;                                                      % Update assignment
        z(:,sum(z)==0) = [];                                                % Remove any empty components
    end
    Z{t} = z;                                                               % Save result        
end