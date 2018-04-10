
% compute steady state and the predictions

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/cleanData/';
listData = dir([dataDir '*.mat']);
dataOut = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';

for ff=1:length(listData)
    
    clear cleanData; clear Axx; clear cfg;
    load([dataDir listData(ff).name]);
    
    % reref to average
    cfg.reref         = 'yes';
    cfg.refchannel    = {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    [cleanData] = ft_preprocessing(cfg,cleanData);
    
    % compute steady state
    cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy
    
    %     % hierarchical structure
    %     % create folder for each subject
    %     newFolder = sprintf('%s', [dataOut listData(ff).name(11:13)]);
    %     if ~exist(newFolder, 'dir')
    %         mkdir(newFolder);
    %     end
    %     allcond = unique(cleanData.trialinfo(:,1));
    %     for cond=1:length(allcond)
    %         cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
    %         Axx = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
    %         save([newFolder '/Axx_' listData(ff).name(11:13) '_c' num2str(cond,'%.2d')],'Axx')
    %     end
    
    % 1 Axx per sbj
    allcond = unique(cleanData.trialinfo(:,1));
    for cond=1:length(allcond)
        cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
        [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
    end
    
    
    % get the condition labels
    testedFreq = [85/8 85/16 85/32]; % in Hz this is the onset of the single stimulus
    onTime = [1 2 4 6 7];
    condNb = 1;
    for testFq=1:length(testedFreq)
        for tt = 1:length(onTime)
            Axx(condNb).condLabel = [num2str(testedFreq(testFq),'%.1f') 'Hz ' num2str(onTime(tt)/8*100) '% DC'];
            condNb = condNb+1;
        end
    end
    
    testedFreq = [85/32 85/16 85/8 85/16 85/32]; % in Hz this is the onset of the single stimulus
    onTime = [1 2 4 6 7];
    for tt = 1:length(testedFreq)
        Axx(condNb).condLabel = ['motion ' num2str(testedFreq(tt),'%.1f') 'Hz ' num2str(onTime(tt)/8*100) '% DC'];
        condNb = condNb+1;
    end
    
    Axx(condNb).condLabel = ['motion ' num2str(85/32,'%.1f') 'Hz 50% DC'];
    Axx(condNb+1).condLabel = ['motion ' num2str(85/16,'%.1f') 'Hz 50% DC'];
    save([dataOut 'Axx_' listData(ff).name(11:13)],'Axx')
    
end

%%%%% had everything in one big struct but not very practical
% for ff=1:length(listData)
%     clear cleanData;
%     load([dataDir listData(ff).name]);
%
%     % compute steady state
%     cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy
%
%     allcond = unique(cleanData.trialinfo(:,1));
%     for cond=1:length(allcond)
%         cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
%         [Subj(ff).Axx(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
%     end
%
%     % call function for prediction save it in the same Subj structure
%     Subj(ff).motionPrediction = predictMotion(Subj(ff).Axx);
%     Subj(ff).simultPrediction = predictSimult(Subj(ff).Axx);
% end
%
% save('Subj','Subj')
%
%
% %%% plot results
% load Subj;
% cfg.layout = 'biosemi128.lay';
% cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
% interactiveSteadyStatePlot(cfg,Subj(1).Axx)
% interactiveSteadyStatePlot(cfg,Subj(1).motionPrediction)
% interactiveSteadyStatePlot(cfg,Subj(1).simultPrediction)
%
% interactiveSteadyStatePlot(cfg,[Subj(1).Axx Subj(2).Axx Subj(3).Axx Subj(4).Axx Subj(5).Axx Subj(6).Axx])
% interactiveSteadyStatePlot(cfg,[Subj(1).motionPrediction Subj(2).motionPrediction])
% interactiveSteadyStatePlot(cfg,[Subj(1).simultPrediction Subj(2).simultPrediction])


%%% average across subj

