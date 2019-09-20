function z = nmi2(x, y)
% Compute nomalized mutual information I(x,y)/sqrt(H(x)*H(y)).
%
% This is second version of nmi, does require the cluster sets to have the same
% number of clusters.
%
% Written by Michael Chen (sth4nth@gmail.com).
% Edited by Jeffrey Chan, 4/2013.

    % resize membership vectors to compact length
    assert(numel(x) == numel(y));
    n = numel(x);
    x = reshape(x,1,n);
    y = reshape(y,1,n);
    
    % relabel vertex->partition mapping to start from 1 (the subtraction of the
    % min part)
    l = min(min(x),min(y));
    x = x-l+1;
    y = y-l+1;
    k = max(max(x),max(y));

    % construct membership matrices, n (number of elements) by k (number of
    % clusters)
    idx = 1:n;
    Mx = sparse(idx,x,1,n,k,n);
    My = sparse(idx,y,1,n,k,n);
    Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
    Hxy = -dot(Pxy,log2(Pxy+eps));

    Px = mean(Mx,1);
    Py = mean(My,1);

    % entropy of Py and Px
    Hx = -dot(Px,log2(Px+eps));
    Hy = -dot(Py,log2(Py+eps));

    % mutual information
    MI = Hx + Hy - Hxy;

    % normalized mutual information
    z = sqrt((MI/Hx)*(MI/Hy));
end