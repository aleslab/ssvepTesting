
% check accuracy perf in Duty Cycle
% for each participant
clearvars
close all

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/originalData/';
behavFiles = dir([dataDir '*.mat']);


keepSBJ = 1:length(behavFiles);



for ss=1:length(keepSBJ)
    ff = keepSBJ(ss);
    clear tableData corResp response corResp matResponse
    load([dataDir behavFiles(ff).name])
    
    % remove any invalid trial from the data
    if length(experimentData) ~= 242 || ss==8 % this sbj has an invalid trial but still 242 trials (restarted many times)
        toRemove = [];
        for ll=1:length(experimentData)
            if experimentData(ll).validTrial == 0 || isempty(experimentData(ll).trialData)
                toRemove = [toRemove; ll];
            end
        end
       keepIndexes = setdiff(1:length(experimentData), toRemove);
       tableData = struct2table(experimentData(keepIndexes));
    else
       tableData = struct2table(experimentData);
    end
    
    %%% analysis per condition - but not many/enough trial per condition 
    %%% for the confusion matrix
    allCond = unique(tableData.condNumber);
    for cc=1:length(allCond)
        indexCond = find(tableData.validTrial & tableData.condNumber == allCond(cc));
        for ic = 1: length(indexCond)
            response(ic,1) =  tableData.trialData(indexCond(ic)).nbDots;
            response(ic,2) = tableData.trialData(indexCond(ic)).response;
            corResp(cc,ic) = tableData.trialData(indexCond(ic)).correct;
        end
        posDots = unique(response(:,1)) + 1;
        for dots = 1: length(posDots)
            trialPerCond = length(find(response(:,1) == dots-1));
            for nbResp = 1: length(posDots)
                matResponse(cc,dots,nbResp) = length(find(response(:,1) == dots-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
            end
        end
    end
%     % plot
%     for cc=1:length(allCond)
%         figure;hold on;
%         imagesc(squeeze(matResponse(cc,:,:)));colorbar;
%     end
    figure;bar(mean(corResp,2)); ylim([0 1]);
    behav(ss,:) = mean(corResp,2);
    
    %%% analysis across conditions
    clear corResp response corResp matResponse
    allCond = unique(tableData.condNumber);
    indexCond = find(tableData.validTrial);
    if ss==8
    for ic = 1: length(indexCond)
        response(ic,1) =  tableData.trialData{indexCond(ic)}.nbDots;
        response(ic,2) = tableData.trialData{indexCond(ic)}.response;
        corResp(ic) = tableData.trialData{indexCond(ic)}.correct;
    end
    else
    for ic = 1: length(indexCond)
        response(ic,1) =  tableData.trialData(indexCond(ic)).nbDots;
        response(ic,2) = tableData.trialData(indexCond(ic)).response;
        corResp(ic) = tableData.trialData(indexCond(ic)).correct;
    end
    end
    posDots = unique(response(:,1)) + 1;
    for nbDim = 1: length(posDots)
        trialPerCond = length(find(response(:,1) == nbDim-1));
        for nbResp = 1: length(posDots)
            matResponse(nbDim,nbResp) = length(find(response(:,1) == nbDim-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
        end
    end
%     % plot
%     figure;imagesc(matResponse);colorbar;
    
end

% each participant overall perf
mean(behav,2)
% average perf per cond
figure;bar(mean(behav)); ylim([0 1]);

