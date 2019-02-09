
% check accuracy perf in Duty Cycle
% for each participant
clearvars
close all

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/originalData/';
behavFiles = dir([dataDir '*.mat']);


keepSBJ = 1:length(behavFiles);
keepSBJ = [1:5 7:11 13:15]; 



for ss=1:length(keepSBJ)
    ff = keepSBJ(ss);
    clear tableData corResp response corResp
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
                matResponse(ss,cc,dots,nbResp) = length(find(response(:,1) == dots-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
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
    clear corResp response corResp
    allCond = unique(tableData.condNumber);
    indexCond = find(tableData.validTrial);
    for ic = 1: length(indexCond)
        response(ic,1) =  tableData.trialData(indexCond(ic)).nbDots;
        response(ic,2) = tableData.trialData(indexCond(ic)).response;
        corResp(ic) = tableData.trialData(indexCond(ic)).correct;
    end
    posDots = unique(response(:,1)) + 1;
    for nbDim = 1: length(posDots)
        trialPerCond = length(find(response(:,1) == nbDim-1));
        for nbResp = 1: length(posDots)
            matResponse2(nbDim,nbResp) = length(find(response(:,1) == nbDim-1 & response(:,2) == nbResp-1)) / trialPerCond *100;
        end
    end
%     % plot
%     figure;imagesc(matResponse);colorbar;
    
end

% each participant overall perf
mean(behav,2)
% average perf per cond
figure;bar(mean(behav)); ylim([0 1]);

% % CI = mean +/- 1.96 * std / sqrt(n)
% ci = 1.96 * std(behav) / sqrt(size(behav,1));

%% plot table behaviour
tabMean = [mean(behav(:,1:5));mean(behav(:,6:10));mean(behav(:,11:15))];
tabMeanMot = repmat(0,3,5);
tabMeanMot(3,1)= mean(behav(:,16));
tabMeanMot(2,2)= mean(behav(:,17));
tabMeanMot(1,3)= mean(behav(:,18));
tabMeanMot(2,4)= mean(behav(:,19));
tabMeanMot(3,5)= mean(behav(:,20));
tabMeanMot(3,3)= mean(behav(:,21));
tabMeanMot(2,3)= mean(behav(:,22));
figure;
subplot(1,2,1)
imagesc(tabMean,[0 1]);colorbar;
subplot(1,2,2)
imagesc(tabMeanMot,[0 1]);colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 1000 1200 500])
saveas(gcf,'detectionPerf.png')



testedFreq = [85/8 85/16 85/32]; % in Hz this is the onset of the single stimulus
onTime = [1 2 4 6 7];
cc=1;
for tt = 1:length(testedFreq)
    for oo=1:length(onTime)
        stimTab(tt,oo) = onTime(oo)/8/testedFreq(tt);
        stimOn(cc) = onTime(oo)/8/testedFreq(tt);
        cc=cc+1;
    end
end
testedFreq = [85/32 85/16 85/8 85/16 85/32 85/32 85/16];
onTime = [1 2 4 6 7 4 4];cc=1;
for tt = 1:length(testedFreq)
    stimOnMov(cc) = onTime(tt)/8/testedFreq(tt);
    cc=cc+1;
end
imagesc(stimOn)

figure;hold on;
for pp=1:3
    scatter(stimOn(5*(pp-1)+1:5*(pp-1)+5),mean(behav(:,5*(pp-1)+1:5*(pp-1)+5)),'filled')
end
scatter(stimOnMov,mean(behav(:,16:22)),'g','filled')
ylabel('accuracy')
xlabel('stim On')
legend('10','5','2','moving','Location','Best')
saveas(gcf,'accByStimOn.png')