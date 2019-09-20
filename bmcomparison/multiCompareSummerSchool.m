function [cmDistances] = multiCompareSummerSchool(csNames)
%
% Compares all the pairwise possibiliteis between among the vsNames
%


sGraphFilename = '/home/jefcha/Research/Data/summerSchool/survey_numeric_combined_edge_list.mat.csv';
% cMeasures = { {'edit', '', ''}, {'', 'jaccard', ''}, {'','rand',''}, {'','vi',''}, {'edgeRecon','',''}, {'pearsonRecon','',''}, {'', 'esVi',''},{'edgeRecon','vi','linear'}};
cMeasures = { {'', 'jaccard', ''}, {'','rand',''}, {'','vi',''}, {'edgeRecon','',''}, {'', 'esVi',''},{'edgeRecon','vi','linear'}};
 
cmDistances = cell(1, size(cMeasures,2));
    
for m = 1 : size(cMeasures,2)
    vMeasure = cMeasures{m}

    mDistances = zeros(size(csNames,2), size(csNames,2));
    
for i = 1 : size(csNames,2)
    sFileName1 = strcat('/home/jefcha/Research/Results/bmcompare/summerSchool/', csNames{i})
    for j = i+1 : size(csNames,2)
        sFileName2 = strcat('/home/jefcha/Research/Results/bmcompare/summerSchool/', csNames{j})
        vDist = compareOverlapBMTemp(sFileName1, sGraphFilename, sFileName2, sGraphFilename, true,true,true,73,73,false,true,true, vMeasure{1}, vMeasure{2}, vMeasure{3},'temp.dis.csv',0.5);
        mDistances(i,j) = vDist(1);
    end 
end

    cmDistances{m} = mDistances;

end % end of for of cMeasures


end