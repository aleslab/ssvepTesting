
% compute steady state

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/MAE/cleanData/';
listData = dir([dataDir '*.mat']);
dataOut = '/Users/marleneponcet/Documents/data/MAE/Axx/';

for ff=3:length(listData)
    
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
% [statusChannel] = ft_steadystateanalysis(cfg, cleanData);
    end
        
    
    % get the condition labels
    spatialFq = [0.125 2]; % standard fq 0.5 /4 or *4
    phase = [10 90];
    dirt = {'left' 'right' 'none'};  % of standard 0=left, 180=right, 99=none (no drift)
    condNb = 1;
    for aa = 1:length(dirt)
        for bb = 1:length(spatialFq)
            for cc=1:length(phase)
                Axx(condNb).condLabel = [dirt{aa} ' sf' num2str(spatialFq(bb)) ' p' num2str(phase(cc))] ;
                condNb = condNb+1;
            end
        end
    end
    save([dataOut 'Axx_' listData(ff).name(1:7)],'Axx','cfg')
    
    
    
    
    test = {'Stat' 'Dyn'};
    stim = [1 2];
    adapt = {'left' 'right' 'none'};
    condNb = 1;
    for aa = 1:length(adapt)
        for ss = 1:length(stim)
            for cc=1:length(test)
                Axx(condNb).condLabel = [num2str(stim(ss)) 'stim' test{cc} 'adap' adapt{aa}] ;
                condNb = condNb+1;
            end
        end
    end
    save([dataOut 'Axx_MAE_S'],'Axx','cfg')

%     phase = [10 90];
%     dirt = {'left' 'right' 'none'};  % of standard 0=left, 180=right, 99=none (no drift)
%     condNb = 1;
%     for aa = 1:length(dirt)
%             for cc=1:length(phase)
%                 Axx(condNb).condLabel = [dirt{aa} ' p' num2str(phase(cc))] ;
%                 condNb = condNb+1;
%             end
%     end
%     save([dataOut 'Axx_MAE_S17'],'Axx','cfg')


%     condLabels = {'rightH90','rightH10','rightH10double','leftH90','leftH170','leftH170double',...
%         'noH90','noH10','noH10double','leftL90','leftL10','leftL10double'};
%     for cc=1:length(condLabels)
%         Axx(cc).condLabel = condLabels{cc} ;
%     end
%     save([dataOut 'AxxMAEf_S01_v3'],'Axx','cfg')
   



    
end



labels = {'2leftD','2leftS','2RightD','2RightS','1leftS','1rightD','2noneS','2noneD'};
for ll=1:length(labels)
    Axx(ll).condLabel = labels{ll};
end

% %%%%%%
% for cond=1:length(allcond)
%     allTrials = find(cleanData.trialinfo(:,1) == allcond(cond));
%     cfg.trials = allTrials(1:floor(end/2));    
%     [AxxFirst(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
%     cfg.trials = allTrials(ceil(end/2):end) ;   
%     [AxxSecond(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
% end
%     
%     save([ 'AxxFirst' ],'AxxFirst','cfg')
%     save([ 'AxxSecond'],'AxxSecond','cfg')
    
      
