function [mAdjMat, vVertRole, mShuffledAdjMat, vShuffledVertRole] = genBlockmodel(sRoleDist, roleNum, graphSize, bShuffle, sGen, varargin)
%
% Generate blockmodels, using a stocastic blockmodel (generative) model.
% It can also generate background null model edges (see Newman's paper_.
%
% @author Jeffrey Chan
% @Date 11/2011
%
% sRoleDist - probability distribution of the sizes of each position/role
% roleNum - number of roles to generate
% graphSize - number of vertices in the graph generated
% bShuffle - whether to shuffle the generated adjacency matri
% sGen - string specifiying whether to generate blockmatrix ('gen' or
% 'nogen')
% mImageGraph - image matrix, each entry being 0-1, specifying the
% probability (i.e., the density of each block) of an edge existing.  The
% mImageGraph must be a square matrix and the number of rows (columns) =
% size(vRoleProb,2).
%
% mAdjMat - (unshuffled) generated adajency matrix.
% vVertRole - mapping of vertices to roles.
% mShuffledAdjMat - shuffled adajency matrix.
% vShuffledVertRole - shuffled vertices to roles mapping.
%

mystream = RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setDefaultStream(mystream);

optargin = size(varargin,2);
stdargin = nargin - optargin;

if stdargin < 4
    disp(' genBlockmodel(vRoleProb, graphSize, bShuffle, sGen, varargin)');
    exit;
end

currIndex = 1;
% file name
bFileOutput = false;
sFilename = '';

% if image matrix missing, we generate it
if optargin > 0
    if strcmp(sGen, 'gen')
        % generate
        % defaults for now
        blockNum =  roleNum;
        probDense = varargin{1};
        meanDense = varargin{2};
        varDense = varargin{3};
        meanSparse = varargin{4};
        varSparse = varargin{5};
        currIndex = currIndex + 5;
        mImageGraph = genImageMatrix(probDense, blockNum, meanDense, varDense, meanSparse, varSparse);
    else
        % read from arguments
        mImageGraph = varargin{1};
        currIndex = currIndex + 1;
    end
    % more command line arguments
    if optargin >= currIndex
        if strcmp(varargin{currIndex}, 'file')
            bFileOutput = true;
            sAdjFilename = varargin{currIndex+1};
            sRoleFilename = varargin{currIndex+2};
            sImageFilename = varargin{currIndex+3};
        end
    end
end

% seed the randomiser
%rand('seed', sum(100*clock));




% generate role sizes
if strcmp(sRoleDist, 'uniform')
    vRoleProb = zeros(1,roleNum) + 1/roleNum;
elseif strcmp(sRoleDist, 'linear')
    vRoleProb = [roleNum:-1:1]/sum([1:1:roleNum]);
elseif strcmp(sRoleDist, 'random')
    % generate a random sizes that sum up to graphsize, then we divide by
    % graphsize
    vRoleProb = (randfixedsum(roleNum, 1, graphSize, 1, graphSize))' / graphSize;
else
    disp('Invalid role distribution');
    disp(' genBlockmodel(sRoleDist, graphSize, sGen, varargin)');
    exit;
end


assert(size(vRoleProb,2) == size(mImageGraph,2));
assert(size(mImageGraph,1) == size(mImageGraph,2));

% edge value
edgeVal = 1;

% generate membership of the vertices
vRoleNum = mnrnd(graphSize, vRoleProb);



% if there are any zero elements, we random sample one of partitions that
% has more than one element and move an element to each zero element
% partition
vZeroSizedRoles = find(vRoleNum == 0);
for zeroi = 1 : size(vZeroSizedRoles,2)
    vNZSizedRoles = find(vRoleNum > 1);
    nzChosen = randsample(vNZSizedRoles, 1);   
    vRoleNum(nzChosen) = vRoleNum(nzChosen) - 1;
    vRoleNum(vZeroSizedRoles(zeroi)) = vRoleNum(vZeroSizedRoles(zeroi)) + 1;
end

% make sure all roles have at least one element
assert(isempty(find(vRoleNum == 0)));
 
% vertex to role mapping
vVertRole = zeros(1, graphSize);
 
vertIndex = 1;
for i = 1 : size(vRoleNum,2)
    % assign vertex, in order, to each role, based on the numbers from
    % vRoleNum
    for j = 1 : vRoleNum(i)
        vVertRole(vertIndex) = i;
        vertIndex = vertIndex + 1;
    end
end





% adjacency matrix generated
mAdjMat = zeros(graphSize, graphSize);
currRowShift = 1;
currColShift = 1;

% from block generate the number of edges (binomial(p,n), p = density of
% block/probability of an edge, n = number of elements in the block)
for r = 1 : size(mImageGraph,1)
    for c = 1 : size(mImageGraph,2)
        edgeNum = binornd(vRoleNum(r) * vRoleNum(c), mImageGraph(r,c));
        vEdges = randsample(vRoleNum(r) * vRoleNum(c), edgeNum);
        for e = 1 : size(vEdges,1)
            row = floor((vEdges(e) - 1) / vRoleNum(c));
            col = floor(mod((vEdges(e)-1), vRoleNum(c)));
            mAdjMat(currRowShift + row, currColShift + col) = edgeVal;
        end
        
        currColShift = currColShift + vRoleNum(c);
    end
    
    currRowShift = currRowShift + vRoleNum(r);
    currColShift = 1;
end



shuffleNum = 50;
% shuffle adjacency matrix (we need to shuffle the actual partitions also)
% shuffle x number of times
vShuffledVertRole = vVertRole;
mShuffledAdjMat = mAdjMat;
if bShuffle
    for s = 1 : shuffleNum
        vShuffle = randperm(graphSize);
        % shuffle adjacency matrix
        mShuffledAdjMat = mShuffledAdjMat(vShuffle, vShuffle);
        % shuffle mapping
        vShuffledVertRole = vShuffledVertRole(vShuffle);
    end
end

% if file output, we assume it is a script and function will be terminated
% after all processing has been finished
if bFileOutput
    dlmwrite(sAdjFilename, mShuffledAdjMat, 'precision', '%d');
    
    % loop through the partitions and write out to file (note it appends)
    for p = 1 : roleNum
        vIndices = find(vShuffledVertRole == p);
        dlmwrite(sRoleFilename, vIndices - 1, '-append'); 
    end
    
    
    dlmwrite(sImageFilename, mImageGraph, 'precision', '%.2f');
       
    exit
end


end % end of function


function [mImageMatrix] = genImageMatrix(probDense, blockNum, meanDense, stdDense, meanSparse, stdSparse)
%
% Generate image matrix.  Each block is either dense or sparse, depending
% on probDense (1 - probDense) respectively.  The density of each block is
% then generated by a gaussian(meanDense, stdDense) for dense (and for
% sparse, we use meanSparse and varSparse).
%
% probDense - probability that a block is dense
% blockNum - number of blocks
% meanDense - 
%

mImageMatrix = zeros(blockNum, blockNum);
mDense = rand(blockNum, blockNum);

for r = 1 : blockNum
    for c = 1 : blockNum
        if mDense(r,c) < probDense
            density = normrnd(meanDense, stdDense);
            if density >= 1.0, density = 1.0; end
            mImageMatrix(r,c) = density;
        else
            density = normrnd(meanSparse, stdSparse);
            if density <= 0.0, density = 0.0; end            
            mImageMatrix(r,c) = density;
        end
    end
end

end % end of function