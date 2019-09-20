clc;

sPosDist = 'uniform';
sRoleDist = 'uniform';
sDegDist = 'poisson';
%posNum = 2;
roleNum = 4;
graphSize = 300;
bShuffle = true;
%vBackgroundProp = [0.1];
degPara = 2;
minDegree = 2;
maxDegree = 10;
backgroundProp = 0;
sBlockmodelType = 'community';
[mAdjMat, vVertRole, mShuffledAdjMat, vShuffledVertRole,mImageGraphPlanted] = genNewmanBlockmodel(sRoleDist, roleNum, graphSize, sDegDist, degPara, minDegree, maxDegree, backgroundProp, sBlockmodelType, bShuffle);
disp('done');
%[cmAdjMat, vVertPos, mVertMembership, cmImageGraph, cmShuffledAdjMat, vShuffledVertPos] = genBlockmodelWithBackground(sPosDist, posNum, graphSize, bShuffle, vBackgroundProp, sBlockmodelType)