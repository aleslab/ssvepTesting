function preProcessing_MAE(sbjNb)
% 12 conditions * 2 trials per condition
% trial: 30s adaptation 5s followed by test 10s adapt 5s test * 8
% adapt frequency = 5 Hz, test frequency = 4.05 Hz (85/21)
% stimulus presented below fixation cross

%%% 1st sbj: nb of frames per ON/OFF is wrong
% might want to change 4Hz to 4.25 Hz 85/20

ff = sbjNb;
clear data, clear cfg, clear cleanData;

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

% path to the files
dataDir = '/Users/marleneponcet/Documents/data/MAE/originalData/';
dataOut = '/Users/marleneponcet/Documents/data/MAE/cleanData/';
eegFiles = dir([dataDir '*.bdf']);
behavFiles = dir([dataDir '*.mat']);

cfg.dataset   =  [dataDir eegFiles(ff).name];
    
% read the behavioural file to get parameters for SSVEP
load([dataDir behavFiles(ff).name])
cfg.trialdef.freqTag = experimentData(1).condInfo.testFreq; % tagging freq (4.05 Hz = 85/21)
cfg.trialdef.cycleLength = 1/cfg.trialdef.freqTag; % duration one cycle (0.2471)
cfg.trialdef.trialLength = cfg.trialdef.cycleLength*experimentData(1).condInfo.testDuration; % duration test stimulus
    
% define trials
cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 112]; % all possible conditions
cfg.trialdef.testTrigger = 150; % trigger sent at the beginning of test
cfg.trialdef.ssvepTagVal = 1;
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'df_MAE';
[cfg] = ft_definetrial(cfg);
    
% pre-processing
cfg.demean        ='no'; % useless with detrend. applies baseline correction to the entire trial by default (or to the specified window in cfg.baselinewindow)
cfg.reref         = 'yes';
cfg.refchannel    = {'A1'}; % A3 = CPz / use Cz = A1
%    cfg.refchannel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
cfg.lpfilter  = 'yes';
cfg.lpfreq = 85; % screen 85Hz ...  49 would really clear noise, not 85
cfg.hpfreq = 1;
cfg.hpfilter = 'no'; % does NOT work anyway
cfg.detrend = 'yes';
[data] = ft_preprocessing(cfg);

    
%%%%%%%%%% Electrode exchanged
for dd=1:length(data.trial)
    tmpChan12 = data.trial{1,dd}(12,:);
    tmpChan29 = data.trial{1,dd}(29,:);
    data.trial{1,dd}(12,:) = tmpChan29;
    data.trial{1,dd}(29,:) = tmpChan12;
end


% resample the data so we have an integer number of samples per cycle
% and define the trials (trl) based on the resampled data
cfg.newFs = 85*6; %Integer number of samples per monitor refresh (~500)
cfg.trialdef.epochLength = 28/85*6; % size of the window for cutting trials (in seconds)
cfg.trialdef.preStimDuration = 14/85*6 ; %21/85*4
data = resample_ssvep(cfg,data);


save(['/Users/marleneponcet/Documents/data/MAE/' eegFiles(ff).name(1:end-4) '_pilote'],'data','cfg')



    %%%%%%%%%%%%%%%%%%%
    % artefact rejection
    % first check for extrem channels/trials
        cfg.layout = 'biosemi128.lay';
        cfg.method = 'distance';
        cfg.neighbours = ft_prepare_neighbours(cfg, data);

        cfg.method = 'summary';
        cfg.keepchannel = 'repair'; % had to modify ft_rejectvisual line 343 so that layout was taken into account
        [data] = ft_rejectvisual(cfg, data); % if I want to change the way how the channels are interpolated then will have to do channel repair separately (will also not have to change the rejectvisual function)

    %     %Eye blink detection and rejection
    %     cfg.artfctdef.eog.trlpadding   = 0.1;
    %     cfg.artfctdef.eog.fltpadding   = 0.1;
    %     cfg.artfctdef.eog.artpadding   = 0.1;
    %     cfg.artfctdef.eog.channel = artefactChannel{ff};  
    %     cfg.artfctdef.eog.interactive = 'yes'; % threshold,
    %     % cfg.artfctdef.eog.bpfreq     = [1 15];
    %     % cfg.artfctdef.eog.bpfiltord  = 4;
    %     [cfg, artifact] = ft_artifact_eog(cfg,data);
    %     cfg.artfctdef.reject = 'complete';
    %     [cleanData] = ft_rejectartifact(cfg, data);
    %
    % %     % Eye movement detection
    % %     for tt=1:length(data.trial)
    % %         diff(tt) = max(abs(data.trial{1,tt}(131,:)-data.trial{1,tt}(132,:))); % EXG3 and 4
    % %     end
    % %     find(diff>80);
    
    % do the cleaning on 16 channels cap + 4 electrodes around the eyes
    cfg.channel =  {'C28','C29','C30','C16','C17','C18','C4','C21','D4','D23','D19','A1','B22','B26','B4','A19','A7','A15','A23','A28','EXG2', 'EXG3','EXG4'}; % {'fp1','fp2','p4','fz','f3','t7','c3','cz','c4','t8','p4','pz','p3','o1','oz','o2'}
    cfg.viewmode = 'vertical';
    cfg = ft_databrowser(cfg,data);
    cfg.artfctdef.reject = 'complete';
    [cleanData] = ft_rejectartifact(cfg, data);
    
    save([dataOut eegFiles(ff).name(1:end-4) '_clean'],'cleanData')
% end

end
%
% % check EXG channels
% trialNb = 605;
% figure;plot(data.trial{1,trialNb}(130,:),'r')
% hold on;plot(data.trial{1,trialNb}(131,:),'g')
% hold on;plot(data.trial{1,trialNb}(132,:),'k')
%
% figure; hold on;
% for tt=1:length(data.trial)
%     plot(data.trial{1,tt}(131,:)-data.trial{1,tt}(132,:))
% end
% for tt=1:length(data.trial)
%     diff(tt) = max(abs(data.trial{1,tt}(131,:)-data.trial{1,tt}(132,:)));
% end

