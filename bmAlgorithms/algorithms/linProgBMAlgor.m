function [mBestImage, mBestMembership, bestObjVal] = linProgBMAlgor(mAdj, k, convEpsilon, runNum, sUpdateFunc)
%
% Input:
% mAdj      - Adjacency matrix (n x n).
% k         - Number of positions to infer.
% convEpsilon   - Epsilon value to denote convergence.
% runNum    - Number of runs to perform.  Each run we start with a random
%             initial mImage.
% sUpdateFunc - The name of the update function to use.
%
% Output:
% mImage    - Image matrix (k x k).
% mMembership   - Membership matrix (n x k).
%
% Computes the hard approximation of A of blockmodel approximate mMem *
% mImage * mMem'.
% Implements our integer/binary k-meansversio.  Has an random mImage initialisation, like
% in k-means/greedy algorithms, and output the best, in terms of objective,
% from all the runs.
%


% vF = linear coefficients of objective.
% mA = coefficients for inequality constraints
% vB = RHS coefficients of inequalit constraints
[vF, mA, vB] = constructLinearProg(mAdj, k);

% run linear program
% vX = solution of linearised problem
vX = bintprog(vF, mA, vB);

% extract out the relevant factors from vX


end % end of function binarykMeansBMAlgor()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

function [vF, mA, vB] = constructLinearProg(mAdj, k)
%
% Construct the linear program.
%
% vF = linear coefficients of objective.
% mA = coefficients for inequality constraints
% vB = RHS coefficients of inequalit constraints
%

% vX, the optimsed variables, should be in the order of C, M then Y, the
% constructed variable for CMC'
% row-wise order

% square adjacency matrix
assert(size(mAdj,1) == size(mAdj,2));

vertNum = size(mAdj,1);
vF = zeros(vertNum * k + k * k + vertNum^2 * k^4, 1);
% all coefficients for C are zero
% assign mAdj_{i,j}




end % end of function




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [objVal] = euclideanObjective(mAdj, mImage, mMembership)
%
% Computes the Euclidean distance approximation objective.
%

%objVal = norm(mAdj - mMembership * mImage * mMembership',2);
mApprox = mAdj - mMembership * mImage * mMembership';
objVal = trace(mApprox * mApprox');

% convert to number if in sparse representation.
if issparse(objVal)
    objVal = full(objVal);
end

end % end of function
