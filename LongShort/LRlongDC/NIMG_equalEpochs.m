
% compute steady state and the predictions with EQUAL number of epochs in
% each condition
% nb of epochs per condition per sbj: 92, 91, 54, 94, 82, 78, 49, 61, 79,
% 65, 82, 90, 83, 47, 71, 58, 74
% 50% rejection = 51... reject S7 S14?? 

clearvars;
addpath C:\Users\Marlene\Documents\git\fieldtrip-aleslab-fork
addpath C:\Users\Marlene\Documents\git\ssvepTesting\svndlCopy
addpath C:\Users\Marlene\Documents\git\ssvepTesting\biosemiUpdated
addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFunctions
ft_defaults

dataDir = 'E:\JUSTIN\data\AM shortlong\E2 LongRangeV2\cleanData\';
listData = dir([dataDir '*.mat']);    
dataOut = 'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E2LRlongDC\AxxEqual\';


for ff=1:length(listData)
    
    clear cleanData; clear Axx; clear cfg;
    load([dataDir listData(ff).name]);
    cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy
    
    % 1 Axx per sbj
    cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    cfg.layout = 'biosemi128.lay';
    allcond = unique(cleanData.trialinfo(:,1));
    
    % get nb of trial per condition
    for cond=1:length(allcond)
        nbTrials(cond) = length(find(cleanData.trialinfo(:,1) == allcond(cond)));
    end
    min(nbTrials)
    for cond=1:length(allcond)
        allTrials = find(cleanData.trialinfo(:,1) == allcond(cond));
        cfg.trials = datasample(allTrials,min(nbTrials),'Replace',false); % select a random subset of trials
        [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
    end
    
    % get the condition labels
    aa = {'Long','Short'}; bb=[0 6];
    for tt=1:2
    Axx(1+bb(tt)).condLabel = [aa{tt} 'Motion'];
    Axx(2+bb(tt)).condLabel = [aa{tt} 'Left'];
    Axx(3+bb(tt)).condLabel = [aa{tt} 'Right'];
    Axx(4+bb(tt)).condLabel = [aa{tt} 'Simult'];
    Axx(5+bb(tt)).condLabel = [aa{tt} 'HalfLeft'];
    Axx(6+bb(tt)).condLabel = [aa{tt} 'HalfRight'];
    end
    
    save([dataOut 'Axx_' listData(ff).name(1:13)],'Axx','cfg')
    
end



