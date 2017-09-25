load('longRange_01__20170727_153334.mat')
% ptbCorgiData = overloadOpenPtbCorgiData('D:\StAndrews\LongRangeSSVEP\longRange_01__20170727_153334.mat')
load('/Users/marleneponcet/Documents/LongRangeSSVEP/longRange_00__20170913_132721.mat')

tableData = struct2table(experimentData);




% analysis per condition
% not enough trials per condition for the confusion matrix
clear corResp response corResp matResponse
allCond = unique(tableData.condNumber);
for cc=1:length(allCond)
    indexCond = find(tableData.validTrial & tableData.condNumber == allCond(cc));
    for ic = 1: length(indexCond)
        response(ic,1) =  tableData.trialData{indexCond(ic)}.dims;
        response(ic,2) = tableData.trialData{indexCond(ic)}.response;
        corResp(cc,ic) = tableData.trialData{indexCond(ic)}.correct;
    end
    posDim = unique(response(:,1)) + 1;
    for nbDim = 1: length(posDim)
        trialPerCond = length(find(response(:,1) == nbDim-1));
        for nbResp = 1: length(posDim)
            matResponse(cc,nbDim,nbResp) = length(find(response(:,1) == nbDim-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
        end
    end
end

for cc=1:length(allCond)
    figure;hold on;
    imagesc(squeeze(matResponse(cc,:,:)));colorbar;
end
figure;bar(mean(corResp(:,:),2))




% analysis across conditions
clear corResp response corResp matResponse
allCond = unique(tableData.condNumber);
indexCond = find(tableData.validTrial);
for ic = 1: length(indexCond)
    response(ic,1) =  tableData.trialData{indexCond(ic)}.dims;
    response(ic,2) = tableData.trialData{indexCond(ic)}.response;
    corResp(ic) = tableData.trialData{indexCond(ic)}.correct;
end

posDim = unique(response(:,1)) + 1;
for nbDim = 1: length(posDim)
    trialPerCond = length(find(response(:,1) == nbDim-1));
    for nbResp = 1: length(posDim)
        matResponse(nbDim,nbResp) = length(find(response(:,1) == nbDim-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
    end
end
        
mean(corResp)
imagesc(matResponse);colorbar;
