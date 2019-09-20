function [mAdjMat, vVertRole, mShuffledAdjMat, vShuffledVertRole,mImageGraphPlanted] = genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, minDegree, maxDegree, backgroundProp, sBlockmodelType, bShuffle)
    %
    % Generate blockmodels, using a stocastic blockmodel (generative) model of
    % Newman.  Undirected (directed is very difficult to generate the image
    % matrix while satisfying in and out deg distributions).
    %
    % @author Jeffrey Chan
    % @Date 11/2012
    %
    % sRoleDist - probability distribution of the sizes of each position/role
    % roleNum - number of roles to generate
    % graphSize - number of vertices in the graph generated
    % sDegDist - probability distribution of the vertice degree (poisson (erdos
    % renyi), powerlaw)
    % degPara - parameter used for the degree distribution
    % backgroundProp - relative proportion of planted vs random
    % sBlockmodelType - structure of planted image matrix (community,
    % core-periphery, hierarchical)
    % bShuffle - whether to shuffle the generated adjacency matrix
    %

    mystream = RandStream('mt19937ar','Seed',sum(100*clock));
    RandStream.setGlobalStream(mystream);

    stdargin = nargin;

    if stdargin < 3
        disp('Usage: genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, backgroundProp, sBlockmodelType, bShuffle)');
    end

    % file name
    bFileOutput = false;

    
    % generate role sizes
    %vVertRole = genVertexPositions(sRoleDist, roleNum, graphSize);
    vVertRole_struct = load('C.mat');
    vVertRole = vVertRole_struct.vVertRole;
    % generate edge degrees
    [vDeg] = genVertexDegree(graphSize, sDegDist, degPara, minDegree, maxDegree);
    
    [vTotalDeg, vProbDeg] = computeTotalDeg(roleNum, vVertRole, vDeg);
    
    % generate the planted blockmodel
    mImageGraphPlanted = genImageGraphPlanted(vTotalDeg, sBlockmodelType, roleNum);
    
    % generte the background, random blockmodel
    mImageGraphBackground = genImageGraphBackground(vTotalDeg, roleNum);
    
    % generate the image matrix
    mImageGraph = (1-backgroundProp) * mImageGraphPlanted + backgroundProp * mImageGraphBackground;

    % generate the graph
    [mAdjMat, ~] = genGraph(vProbDeg, vVertRole, mImageGraph);
    
    
    

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [vVertRole] = genVertexPositions(sRoleDist, roleNum, graphSize)
    %
    % Generate the position assignments for each vertex.
    % Vertex sizes are random generated.
    %

    if strcmp(sRoleDist, 'uniform')
        vRoleProb = zeros(1,roleNum) + 1/roleNum;
    elseif strcmp(sRoleDist, 'linear')
        vRoleProb = [roleNum:-1:1]/sum([1:1:roleNum]);
    elseif strcmp(sRoleDist, 'random')
        % generate a random sizes that sum up to graphsize, then we divide by
        % graphsize
        vRoleProb = (randfixedsum(roleNum, 1, graphSize, 1, graphSize))' / graphSize;
    else
        disp('Usage: genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, backgroundProp, sBlockmodelType, bShuffle)');
        err = MException('genVertexPositions:InvalidRoleDist', 'Invalid role/position distribution');
        throw(err);        
    end
    
    
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
    assert(isempty(find(vRoleNum == 0,1)));
 
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
    
    
end % end of function()


function [vDeg] = genVertexDegree(graphSize, sDegDist, degPara, minDegree, maxDegree)
    %
    % Generate the vertex degree of each vertex.
    %

    switch sDegDist
        case 'powerlaw'
            % generate expected degrees for the vertices
            % P(y) = (y+b)^(-n)
    
            % in degree
            vX = rand(1, graphSize);
            vDeg = vX.^(1/1-degPara);
        case 'poisson'
            vDeg = poissrnd(degPara, 1, graphSize);
        otherwise
            disp('Usage: genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, backgroundProp, sBlockmodelType, bShuffle)');  
            
            err = MException('genVertexDegree:InvalidDegDist', 'Invalid degree distribution');
            throw(err);
    end
    
    % the following could introduce bias into the distributions
    % at least minDegree
    vDeg(vDeg < minDegree) = minDegree;
    % restrict all expected degrees to maxDegree
    vDeg(vDeg > maxDegree) = maxDegree;

end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555

function [vTotalDeg, vProbDeg] = computeTotalDeg(roleNum, vVertRole, vDeg)

    % compute the total in and out degree for each position
    vTotalDeg = zeros(1, roleNum);

    vertNum = size(vDeg,2);
    
    for v = 1 : vertNum
       vTotalDeg(vVertRole(v)) = vTotalDeg(vVertRole(v)) + vDeg(v);
    end
    
    % compute the deg probability 
    vProbDeg = zeros(1, vertNum);
    
    for v = 1 : vertNum
        vProbDeg(v) = vDeg(v) / vTotalDeg(vVertRole(v));
    end
end % end of function


function mImageGraphPlanted = genImageGraphPlanted(vTotalDeg, sBlockmodelType, roleNum)
    %
    % Generate the planted image graph structure.
    %
    

    
    % generate the planted image graph
    mImageGraphPlanted = zeros(roleNum, roleNum);
    switch sBlockmodelType
        case 'community'
            for r = 1 : roleNum
                % we take the average, since both total degrees must be the same
                mImageGraphPlanted(r,r) = vTotalDeg(r);
            end
        case 'coreperiphery'
            %
            % Might throw exception in randfixedsum (line 46) for certain degree distributions
            %
            
            % get the maximum in/out degree as the core
            [vSortedDeg, vSortedIndex] = sort(vTotalDeg, 'descend');
            
            vRemainigRowDeg = vSortedDeg
            % generate the row sums
            for r = 1 : roleNum
                r
                vRemainigRowDeg
                vRowDeg = randfixedsum(roleNum-r+1, 1, vRemainigRowDeg(r), 0, vRemainigRowDeg(r))';
                vSortedRowDeg = sort(vRowDeg, 'descend');               
                for c = r : roleNum
                    mImageGraphPlanted(r,c) = vSortedRowDeg(c-r+1);
                    mImageGraphPlanted(c,r) = vSortedRowDeg(c-r+1);
                    if c == r
                        vRemainigRowDeg(r) = vRemainigRowDeg(r) - mImageGraphPlanted(r,r);
                    else
                        vRemainigRowDeg(r) = vRemainigRowDeg(r) - mImageGraphPlanted(r,c);
                        vRemainigRowDeg(c) = vRemainigRowDeg(c) - mImageGraphPlanted(r,c);
                    end
                end
                
%                 vRemainigRowDeg = vRemainigRowDeg - sum(mImageGraphPlanted, 1);   
            end 
            
            % revert ordering
            mImageGraphPlanted(vSortedIndex, vSortedIndex) = mImageGraphPlanted;
        case 'hierarchical'
            % get the two smallest total degrees
            [vSortedDeg, vSortedIndex] = sort(vTotalDeg, 'ascend');
            
            
            % assign the first and last row using the first two degrees
            mImageGraphPlanted(1,2) = vSortedDeg(1);
            mImageGraphPlanted(2,1) = vSortedDeg(1);
%             mImageGraphPlanted(roleNum, roleNum-1) = vSortedDeg(2);
%             mImageGraphPlanted(roleNum-1, roleNum) = vSortedDeg(2);
            
            % rest
            vRemainingDeg = vSortedDeg(3:end);
            vOrder = randperm(roleNum-2);
            %vRemainingDeg = vRemainingDeg(randperm(roleNum-2));
            for r = 1 : roleNum - 2
                assert(mImageGraphPlanted(r+1, r) > 0);
                mImageGraphPlanted(r+1, r+2) =  vSortedDeg(r+1) - mImageGraphPlanted(r+1, r);
                mImageGraphPlanted(r+2, r+1) = mImageGraphPlanted(r+1, r+2);
            end
            
            % remaining degrees in a diagonal element
            mImageGraphPlanted(roleNum, roleNum) = vSortedDeg(roleNum) - mImageGraphPlanted(roleNum-1, roleNum);
            
            % revert ordering
            mImageGraphPlanted(vSortedIndex, vSortedIndex) = mImageGraphPlanted;
        otherwise
            disp('Usage: genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, backgroundProp, sBlockmodelType, bShuffle)');  
            
            err = MException('genImageGraphPlanted:InvalidBlockType', 'Invalid blockmodel type specified.');
            throw(err);
    end
    

end % end of function


function mImageGraphBackground = genImageGraphBackground(vTotalDeg, roleNum)
    %
    % Compute the background image matrix.
    %
    
    edgeNum = sum(vTotalDeg);
    
    mImageGraphBackground = zeros(roleNum, roleNum);
    
    for r = 1 : roleNum
        for c = 1 : roleNum
            mImageGraphBackground(r,c) = vTotalDeg(r) * vTotalDeg(c) / (2 * edgeNum);
        end
    end
    
    
end % end of function()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


function [mAdjMat, mNormalisedImageGraph] = genGraph(vProbDeg, vVertRole, mImageGraph)
    %
    % Generate the graph and also return the normalised image graph.
    %
    
    posNum = size(mImageGraph,1);
    vertNum = size(vVertRole,2);
    % find the vertices in each position
    cVertProb = cell(posNum, 1);
    % map of each position and the vertices it contains
    cPosVert = cell(posNum, 1);
    
    mNormalisedImageGraph = zeros(posNum, posNum);
    
    for p = 1 : posNum
        [vI, vP] = find(vVertRole == p);
        cVertProb{p} = vProbDeg(vP);
        cPosVert{p} = vP;
    end
    
    mAdjMat = zeros(vertNum, vertNum);
    % generate the number of edges for each block, then allocate the incident
    % vertices 
    % note this generates symmetric graphs
    for r = 1 : posNum
        for c = r : posNum
            edgeNum = mImageGraph(r,c);
            % sample the edges, and if an edge already exist, go again
            % determine size of block
            blockSize = size(cVertProb{r},2) * size(cVertProb{c},2);
            
            if edgeNum > blockSize
                edgeNum = blockSize;
            end
            
            assert(blockSize > 0);
            mNormalisedImageGraph(r,c) = edgeNum / blockSize;
%             assert(edgeNum <= blockSize);
            
            vSrcVerts = cPosVert{r};
            vTarVerts = cPosVert{c};

            % see if we generate edges or non-edges
            if edgeNum <= blockSize / 2
                vCdfRowPos = cumsum(cVertProb{r});
                vCdfColPos = cumsum(cVertProb{c});
               
                genEdgeNum = 0;
                while genEdgeNum < edgeNum
                    rowNum = rand;
                    colNum = rand;
                    srcVert = vSrcVerts(sum(vCdfRowPos < rowNum)+1);
                    tarVert = vTarVerts(sum(vCdfColPos < colNum)+1);
                    % only add edge if it doesn't exist
                    if mAdjMat(srcVert, tarVert) == 0
                        mAdjMat(srcVert, tarVert) = 1;
                        % symmetric
                        mAdjMat(tarVert, srcVert) = 1;
                        genEdgeNum = genEdgeNum + 1;
                    end
               end
            else
                % assign block to dense first
                mAdjMat(vSrcVerts, vTarVerts) = 1;
                % symmetric
                mAdjMat(vTarVerts, vSrcVerts) = 1;
                % invert vertex probabilites
                if length(cVertProb{r}) > 1
                    vCdfRowPos = cumsum(1-cVertProb{r});
                else
                    vCdfRowPos = cumsum(cVertProb{r});
                end
                if length(cVertProb{c}) > 1
                    vCdfColPos = cumsum(1-cVertProb{c});
                else
                    vCdfColPos = cumsum(cVertProb{c});
                end
               
                genEdgeNum = 0;
                while genEdgeNum < blockSize - edgeNum
                    rowNum = rand;
                    colNum = rand;
                    srcVert = vSrcVerts(sum(vCdfRowPos < rowNum)+1);
                    tarVert = vTarVerts(sum(vCdfColPos < colNum)+1);
                    % only add edge if it doesn't exist
                    if mAdjMat(srcVert, tarVert) == 1
                        mAdjMat(srcVert, tarVert) = 0;
                        % symmetric
                        mAdjMat(tarVert, srcVert) = 0;
                        genEdgeNum = genEdgeNum + 1;
                    end
               end                
            end
                        
        end
    end
end % end of function()