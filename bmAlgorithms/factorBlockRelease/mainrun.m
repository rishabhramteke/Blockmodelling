clc;clear;

%% input 
A = csvread('baboons.el.csv')
GT_ = csvread('baboons.part.csv');
%GT = zeros(14,1);
%baboons
for i = 1:8
    GT(GT_(1,i),1) = 1;     
end
for i = 1:6
    GT(GT_(2,i),1) = 2;
end
%polbooks 105x105
%{
GT = zeros(105,1);
for i = 1:49
    GT(GT_(1,i),1) = 1;     
end
for i = 1:43
    GT(GT_(2,i),1) = 2;
end
for i = 1:13
    GT(GT_(3,i),1) = 3;
end
%}


%polblog dataset
%{
GT = zeros(1490,1);
for i = 1:758
    GT(GT_(1,i),1) = 1;     
end
for i = 1:732
    GT(GT_(2,i),1) = 2;
end
%}

%check = load('adjnoun.mat');
%baboons is 14x14
%football datset
%{
GT = zeros(115,1);
for i = 1:9
    GT(GT_(1,i),1) = 1;     
end
for i = 1:8
    GT(GT_(2,i),1) = 2;
end
for i = 1:11
    GT(GT_(3,i),1) = 3;
end
for i = 1:12
    GT(GT_(4,i),1) = 4;
end
for i = 1:10
    GT(GT_(5,i),1) = 5;
end
for i = 1:5
    GT(GT_(6,i),1) = 6;
end
for i = 1:13
    GT(GT_(7,i),1) = 7;
end
for i = 1:8
    GT(GT_(8,i),1) = 8;
end
for i = 1:10
    GT(GT_(9,i),1) = 9;
end
for i = 1:12
    GT(GT_(10,i),1) = 10;
end
for i = 1:7
    GT(GT_(11,i),1) = 11;
end
for i = 1:10
    GT(GT_(12,i),1) = 12;
end
%}


%edgeList = csvread('ptbr_edges.csv', 1, 1);
%adj = sparse(edgeList(1, :), edgeList(2, :), 1, 1978, 1978);
%mAdj=full(adj); % n = 1978

%Extracting edges from gml file graph
%{
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
      if length(nums)
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

%}

adj = sparse(A(:, 1), A(:, 2), 1, 14, 14);
mAdj=full(adj);
k = 2;
convEpsilon = 1e-8 ;
runNum = 10;
sUpdateFunc = 'projGradDescM' ;
sDistanceFunc = 'bmEuclideanAdj';
sUpdateMembershipFunc = 'hardCIncr' ;
sImageInit = 'randomInit' ;
sMemInit = 'randomInit' ;


%% output
[mBestImage, mBestMembership, bestObjVal, totalIterNum] = binaryMembershipBMAlgor(mAdj, k, convEpsilon, runNum, sUpdateFunc, sDistanceFunc, sUpdateMembershipFunc, sImageInit, sMemInit);
C = mBestMembership;
M = mBestImage;
error = norm(mAdj -  C * M * C','fro')
%simi = C' * C * M * C' * C ;
Z = zeros(14,1);
for i = 1:14
    for j = 1:k
        if max(C(i,:)) == C(i,j)
            Z(i,1) = j;
        end
    end
        
end

%for i = 1:70
%    Z(i,1) = GT(i,1);
%end


NMI_1 = nmi2(Z,GT)