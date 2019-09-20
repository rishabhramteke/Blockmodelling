clc;clear;

%% input 
%A = csvread('ptbr_edges.csv', 1, 1)
%edgeList = csvread('ptbr_edges.csv', 1, 1);
%adj = sparse(edgeList(1, :), edgeList(2, :), 1, 1978, 1978);
%mAdj=full(adj); % n = 1978

%Extracting edges from gml file graph
fileName = 'karate.gml';
inputfile = fopen(fileName);
A=[];
l=0;
k=1;
while 1
      % Get a line from the input file
      tline = fgetl(inputfile);
      % Quit if end of file
      if ~ischar(tline)
          break
      end
      nums = regexp(tline,'\d+','match');
      if ~isempty(nums)
          if l==1
              l=0;
              A(k,2)=str2num(nums{1});  
              k=k+1;
              continue;
          end
          A(k,1)=str2num(nums{1});
          l=1;
      else
          l=0;
          continue;
      end
end

adj = sparse(A(:, 1), A(:, 2), 1, 34, 34);
mAdj=full(adj);


k = 2;
convEpsilon = 1e-8 ;
runNum = 10;
sUpdateFunc = 'projGradDescM' ;
sDistanceFunc = 'euclidean';
sUpdateMembershipFunc = 'hardCIncr' ;
sImageInit = 'randomInit' ;
sMemInit = 'randomInit' ;
fImageDistanceFunc = 'euclidean';
cfValidMeasFunc = 'euclidean';
bDiscretiseMembership = false;
bColNormalise = false;

%% output
[mBestImage, mBestMembership, bestObjVal, cvObjVal, cvImageDis, cvKKTResidual, ccvGroundComparison, ccmImage, ccmMembership, totalIterNum] = binaryMembershipBMAlgor(mAdj, k, convEpsilon, runNum, sUpdateFunc, sDistanceFunc, fImageDistanceFunc, sUpdateMembershipFunc, sImageInit, sMemInit, cfValidMeasFunc, bDiscretiseMembership, bColNormalise)