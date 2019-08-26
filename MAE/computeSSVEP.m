
% compute steady state

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/MAE/cleanData/';
listData = dir([dataDir '*.mat']);
dataOut = '/Users/marleneponcet/Documents/data/MAE/Axx/';

for ff=12:length(listData)
    
    clear cleanData; clear Axx; clear cfg;
    load([dataDir listData(ff).name]);
    
    % compute steady state
    cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy
    
    % 1 Axx per sbj
    cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    cfg.layout = 'biosemi128.lay';
    allcond = unique(cleanData.trialinfo(:,1));
        
    for cond=1:length(allcond)
        cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
        [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
    end
        
    
    % get the condition labels
    spatialFq = [0.125 2]; % standard fq 0.5 /4 or *4
    phase = [5 90];
    dir = {'left' 'right' 'none'};  % of standard 0=left, 180=right, 99=none (no drift)
    condNb = 1;
    for aa = 1:length(dir)
        for bb = 1:length(spatialFq)
            for cc=1:length(phase)
                Axx(condNb).condLabel = [dir{aa} ' sf' num2str(spatialFq(bb)) ' p' num2str(phase(cc))] ;
                condNb = condNb+1;
            end
        end
    end
    
    save([dataOut 'Axx_' listData(ff).name(11:21)],'Axx')
    
end


