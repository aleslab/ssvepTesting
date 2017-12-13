
% check accuracy perf in LongRangeV2

dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/originalData/';
behavFiles = dir([dataDir '*.mat']);

for ff=1:length(behavFiles)
    load([dataDir behavFiles(ff).name])


load('longRange_01__20170727_153334.mat')


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
