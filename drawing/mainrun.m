clc;clear;

%% input 
BM = load('Polbooks.mat');
mAdj = BM.C;
mMembership = BM.C;
vVertLineLocs = [];
vHorLineLocs = [];
mPosColour = [0 1 0; 0 0 0];
%% output
%drawColouredBlockmodel(mAdj, mMembership, mPosColour, vVertLineLocs, vHorLineLocs)
drawBlockmodel(mAdj, vVertLineLocs, vHorLineLocs)
