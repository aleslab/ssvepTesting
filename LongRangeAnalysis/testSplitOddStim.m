
load cleanData_S02.mat
load('/Users/marleneponcet/Documents/LongRangeSSVEP/longRange_00__20170913_132721.mat')

% we need to know which cycle (1-5) contains a dim

% there is a total of 62 refresh for a 12.4 s trial
% a 2 s epoch is 10 refresh
% do not include begining and end (just keep the 10 s)
% so we are interested in the 7th frame up to 56 (included)
limits = 7:10:57;
for cc=1:length(limits)-1
    cycle(cc,:) = limits(cc):limits(cc+1)-1;
end


% find out the nb and when the dims are
% and create a new column that determines if the cycle is dim (0) or not
% (1)
for trial=1:length(cleanData.trial)
    originalTrial = cleanData.trialinfo(trial,2);
    allDims = experimentData(originalTrial).trialData.stimDim;
    currentCycle = cleanData.trialinfo(trial,4);
    if find(ismember(allDims,cycle(currentCycle,:))==1)
        cleanData.trialinfo(trial,5) = 0;
    else
        cleanData.trialinfo(trial,5) = 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
allcond = unique(cleanData.trialinfo(:,1));

%%%%%%%%%%%%%%%%%%%%%
%%%%% 'NORMAL' (not sweep)
%%%%%%%%%%%%%%%%%%%%%
for cond=1:8
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond) & cleanData.trialinfo(:,5) == 1);
    [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData);
end
for cond=9:16
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond-8) & cleanData.trialinfo(:,5) == 0);
    [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData);
end
save('testSplitOdd', 'Axx')




%%%%%%%%%%%%%%%%
addpath /Users/marleneponcet/Documents/Git/fieldtrip
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
interactiveSteadyStatePlot(cfg,Axx)
