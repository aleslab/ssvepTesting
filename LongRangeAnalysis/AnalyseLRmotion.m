clear all;

addpath /Users/marleneponcet/Documents/Git/fieldtrip
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults


% path to the files
eegFiles = {'/Users/marleneponcet/Documents/LongRangeSSVEP/LongRangeS01.bdf'
    '/Users/marleneponcet/Documents/LongRangeSSVEP/mp-lrmot-20170913-1.bdf'};

behavFiles = {'/Users/marleneponcet/Documents/LongRangeSSVEP/longRange_01__20170727_153334.mat'
    '/Users/marleneponcet/Documents/LongRangeSSVEP/longRange_00__20170913_132721.mat'};

% channel used for artefact rejection
artefactChannel = {'C17', 'EXG4'};

for ff=1:length(eegFiles)
    clear data, clear cfg, clear cleanData;
    cfg.dataset   =  eegFiles{ff};
    
    % read the behavioural file
    load(behavFiles{ff})
    tableData = struct2table(experimentData);
    
    % save some of the data from the behavioural file into the cfg
    cfg.trialdef.nbTotalCycles = experimentData(1).trialData.nbTotalCycles;
    cfg.trialdef.preStimDuration = experimentData(1).trialData.preStimDuration;
    cfg.trialdef.trialLength = experimentData(1).trialData.trialDuration;
    cfg.trialdef.freqTag = experimentData(1).condInfo.stimTagFreq;
    cfg.trialdef.cycleLength = 1/cfg.trialdef.freqTag;
    
    % define trials
    cfg.trialdef.bitmask = 2^9-1; %Values to keep.
    cfg.trialdef.condRange = [101 165]; % all possible conditions
    cfg.trialdef.ssvepTagVal = 1;
    cfg.layout = 'biosemi128.lay';
    cfg.trialfun = 'LR_fullTrial';
    [cfg] = ft_definetrial(cfg);
    
    % pre-processing
    cfg.lpfilter  = 'yes';
    cfg.lpfreq        = 48; % < screen 85Hz
    cfg.demean        ='yes';
    cfg.reref         = 'yes';
    cfg.refchannel    = {'A1'}; % A3 = CPz / use Cz = A1
    [data] = ft_preprocessing(cfg);
    
    % resample the data so we have an integer number of samples per cycle
    % and define the trials (trl) based on the resampled data
    cfg.newFs = 85*6; %Integer number of samples per monitor refresh (~500)
    cfg.trialdef.epochLength = 2; % size of the window for cutting trials (in seconds)
    data = resample_ssvep(cfg,data);
    
    %%%%%%%%%%%%%%%%%%%
    % artefact rejection
    % first check for extrem channels/trials
    cfg.layout = 'biosemi128.lay';
    cfg.method = 'distance';
    cfg.neighbours = ft_prepare_neighbours(cfg, data);
    
    cfg.method = 'summary';
    cfg.keepchannel = 'repair'; % had to modify ft_rejectvisual line 343 so that layout was taken into account
    [data] = ft_rejectvisual_modif(cfg, data); % if I want to change the way how the channels are interpolated then will have to do channel repair separately (will also not have to change the rejectvisual function)
    
    %Eye blink detection and rejection
    cfg.artfctdef.eog.trlpadding   = 0.1;
    cfg.artfctdef.eog.fltpadding   = 0.1;
    cfg.artfctdef.eog.artpadding   = 0.1;
    cfg.artfctdef.eog.channel = artefactChannel{ff};
    cfg.artfctdef.eog.interactive = 'yes'; % threshold, S1=3, S2=
    % cfg.artfctdef.eog.bpfreq     = [1 15];
    % cfg.artfctdef.eog.bpfiltord  = 4;
    [cfg, artifact] = ft_artifact_eog(cfg,data);
    
    cfg.artfctdef.reject = 'complete';
    [cleanData] = ft_rejectartifact(cfg, data);
    
%     % optional: last look
%     cfg.viewmode = 'vertical';
%     cfg = ft_databrowser(cfg,cleanData);
%     cfg.artfctdef.reject  = 'complete';
%     cleanData = ft_rejectartifact(cfg,cleanData);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cfg.layout = 'biosemi128.lay';
    cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
    allcond = unique(cleanData.trialinfo(:,1));
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%%% 'NORMAL' (not sweep)
    %%%%%%%%%%%%%%%%%%%%%
    for cond=1:8
        cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
        [Subj(ff).Axx(cond)] = ft_steadystateanalysis(cfg, cleanData);
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%%% SWEEP
    %%%%%%%%%%%%%%%%%%%%%
    % change in position every 2 sec
    % use column 4 of trialinfo to know at which position the stimulus is
    % the cycles (5 per 12s trial) and the positions are the same
    % between 15 to 19 epochs per condition...
    for cycle=1:5
        cfg.trials = find(cleanData.trialinfo(:,1) == 109 & cleanData.trialinfo(:,4) == cycle);
        [Subj(ff).Sweep(cycle)] = ft_steadystateanalysis(cfg, cleanData);
    end
    
end
save('Subj','Subj')


%%% plot results
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
interactiveSteadyStatePlot(cfg,[Subj(1).Axx Subj(2).Axx])
interactiveSteadyStatePlot(cfg,[Subj(1).Sweep Subj(2).Sweep])








%%% Cleaning
% S01: channels A13 D7 C7  replaced, trial number 280 removed, 
% 269 trials with blinks removed
% S02: channeks D1, D2, D10 replaced, trials 3, 88, 227, 363, 453 removed, 
% 33 trials with blinks removed

