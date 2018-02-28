
% compute steady state and the predictions

clear vars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/cleanData/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    
dataOut = '/Users/marleneponcet/Documents/dataLongRangeV2/AxxEpoch/';

for ff=13:length(listData)
    fprintf('Processing sbj %s \n',num2str(ff))
    clear cleanData; 
    load([dataDir listData(ff).name]);
    
    % compute steady state
    cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy

    allcond = unique(cleanData.trialinfo(:,1));  
    for cond=1:length(allcond)
        clear Axx;
        % 1st epoch
        cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond) & cleanData.trialinfo(:,4) == 1);
        Axx = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced  
        save([dataOut '/Axx_' listData(ff).name(11:13) '_c' num2str(cond,'%.2d') '_e1'],'Axx')
        % last epoch
        clear Axx;
        cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond) & cleanData.trialinfo(:,4) == 5);
        Axx = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced  
        save([dataOut '/Axx_' listData(ff).name(11:13) '_c' num2str(cond,'%.2d') '_e5'],'Axx')
    end    
     
end




