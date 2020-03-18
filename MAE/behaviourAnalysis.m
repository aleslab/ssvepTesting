
filesDir = '/Users/marleneponcet/Documents/data/MAEv3/original/';
fullList = {dir([filesDir 'MAEv3_*.mat']), dir([filesDir 'MAEv3b_*.mat'])};


% responses:
% 1=left
% 2=right
% 3=no motion


% will need to reject the first answers: will be the same as the previous
% response (= during adaptation)

% 76 cycles during the test, 1 cycle is around 117.6 ms (1/8.5)
% remove the 1st 2 seconds (=17 cycles)? 
% minTime = 18;
minTime = 5;
tabResults = [];
sbjOrderE1 = [105 101 109 102 103 107 104 106 108 110 112 111 113 115 114 116 118 117 120 119];
sbjOrderE2 = [102 104 101 110 106 103 108 107 105 109 111 114 112 115 113 117 116 118 119 120];

for expt = 1:length(fullList)
    filesList = fullList{expt};
    if mod(length(filesList),5)>0
        disp('incorrect number of files !!!!!!!')
    end
    for subj=1:length(filesList)/5
        tabSbj = []; tabTmp=[];
        for ff=1:5
            load([filesDir filesList(ff+5*(subj-1)).name])
            clear condition left right none;
            for trial = 1:length(experimentData)
                condition = experimentData(trial).trialData.trigger;
                left = length(find(experimentData(trial).trialData.allResp(minTime:end) == 1));
                right = length(find(experimentData(trial).trialData.allResp(minTime:end) == 2));
                none = length(find(experimentData(trial).trialData.allResp(minTime:end) == 3));
                total = left + right + none;
                pLeft = left/total; pRight = right/total; pNone = none/total;
                tabTmp = [tabTmp ; table(condition, pLeft, pRight, pNone)];
            end
        end
        % get the average for that sbj
        clear sbjMeanLeft sbjMeanRight sbjMeanNone;
        allCond = unique(tabTmp.condition);
        for iCond=1:length(allCond)
            sbjMeanLeft(iCond) = mean(tabTmp.pLeft(find(tabTmp.condition==allCond(iCond))));
            sbjMeanRight(iCond) = mean(tabTmp.pRight(find(tabTmp.condition==allCond(iCond))));
            sbjMeanNone(iCond) = mean(tabTmp.pNone(find(tabTmp.condition==allCond(iCond))));
        end
        if expt == 1
        tabSbj = table(repmat(subj,length(allCond),1),allCond,sbjMeanLeft',sbjMeanRight',sbjMeanNone');
        else
        tabSbj = table(repmat(subj,length(allCond),1),allCond+30,sbjMeanLeft',sbjMeanRight',sbjMeanNone');
        end
        figure;imagesc([sbjMeanLeft',sbjMeanRight',sbjMeanNone']);
        xlabel({'response'});xticks([1:3]);xticklabels({'left','right','no motion'})
        if expt==1
            ylabel({'condition'});yticklabels({'noG1-10','noG2h-10','noG2h-90','leftG1-10','leftG2h-10','leftG2h-90','rightG1-10','rightG2h-10','rightG2h-90'})
            saveas(gcf,['figures' filesep 'Behav_S' num2str(sbjOrderE1(subj)) 'e' num2str(expt) '.png'],'png')
        else
            ylabel({'condition'});yticklabels({'noG1-90','noG2l-10','noG2l-90','leftG1-90','leftG2l-10','leftG2l-90','rightG1-90','rightG2l-10','rightG2l-90'})
            saveas(gcf,['figures' filesep 'Behav_S' num2str(sbjOrderE2(subj)) 'e' num2str(expt) '.png'],'png')
        end
        tabSbj.Properties.VariableNames = {'sbjNb','condition', 'left', 'right', 'none'};
        tabResults = [tabResults; tabSbj];
    end
end

allCond = unique(tabResults.condition);
for iCond=1:length(allCond)
    meanLeft(iCond) = nanmean(tabResults.left(find(tabResults.condition==allCond(iCond))));
    meanRight(iCond) = nanmean(tabResults.right(find(tabResults.condition==allCond(iCond))));
    meanNone(iCond) = nanmean(tabResults.none(find(tabResults.condition==allCond(iCond))));
end
figure;
imagesc([meanLeft' meanRight' meanNone'])
ylabel({'condition'});yticks([1:18]);
yticklabels({'noG1-10','noG2h-10','noG2h-90','leftG1-10','leftG2h-10','leftG2h-90','rightG1-10','rightG2h-10','rightG2h-90',...
    'noG1-90','noG2l-10','noG2l-90','leftG1-90','leftG2l-10','leftG2l-90','rightG1-90','rightG2l-10','rightG2l-90'})

% re-order for  better understanding
orderCond = [1 4 7 2 5 8 11 14 17 10 13 16 3 6 9 12 15 18];
figure;
imagesc([meanLeft(orderCond)' meanRight(orderCond)' meanNone(orderCond)'])
ylabel({'condition'});yticks([1:18]);
yticklabels({'noG1-10','leftG1-10','rightG1-10','noG2h-10','leftG2h-10','rightG2h-10','noG2l-10','leftG2l-10','rightG2l-10'...
    'noG1-90','leftG1-90','rightG1-90','noG2h-90','leftG2h-90','rightG2h-90','noG2l-90','leftG2l-90','rightG2l-90'});
xlabel({'response'});xticks([1:3]);xticklabels({'left','right','no motion'})
saveas(gcf,['figures' filesep 'Behav_Average.png'],'png')

