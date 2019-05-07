function preProcessing_seqContext_lastStim(sbjNb)
%%% read and extract epochs for seqContext experiment
% 5 conditions * 17 blocks * 2 repetitions per block
% each sequence is repeated 10 times in one trial condition
% fq presentation 85/24 at 50% duty-cycle
% 4 stim followed by 2 stim duration blanks = epochs of 1.7 sec(6*24/85)


%%%%%% ATTENTION
% will probably need to do something different for the epochs in the random
% condition... how to "cut" those???

% also maybe I should only take the last stimulus presention (not the
% entire sequence) = make a different resample_ssvep
% run a sbj with 2 freq?

ff = sbjNb;

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% path to the files
dataDir = '/Users/marleneponcet/Documents/data/seqContext/originalData/';
dataOut = '/Users/marleneponcet/Documents/data/seqContext/cleanData/';
eegFiles = dir([dataDir '*.bdf']);
behavFiles = dir([dataDir '*.mat']);

clear data, clear cfg, clear cleanData;
cfg.dataset   =  [dataDir eegFiles(ff).name];

% read the behavioural file
load([dataDir behavFiles(ff).name])

for trial=1:length(experimentData)
    % oups: have to modify exp program that should not add 0 0 if > 60 stim
    seqStim(trial,:) = experimentData(trial).trialData.fullSeq(1:60);
end

%%%%%%%%%%%%%%%%%%%
%%% define trials (10 repeats of the same sequence)
% ONLY TAKE 2 cycles trial for ERP at the last stim position
%Silly hard coded numbers required for df_seqContext
cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 165]; % all possible conditions
cfg.trialdef.ssvepTagVal = 1;
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'df_seqContext_lastStim';
cfg.stimSeq = seqStim;
[cfg] = ft_definetrial(cfg);

%%% pre-processing
cfg.demean        ='no'; % useless with detrend. applies baseline correction to the entire trial by default (or to the specified window in cfg.baselinewindow)
cfg.reref         = 'yes';
cfg.refchannel    = {'A1'}; % A3 = CPz / use Cz = A1
%    cfg.refchannel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
cfg.lpfilter  = 'yes';
cfg.lpfreq = 85; % screen 85Hz 
cfg.hpfreq = 1;
cfg.hpfilter = 'no'; % does NOT work anyway
cfg.detrend = 'yes';
[data] = ft_preprocessing(cfg);

% exchange A12 and A29 (A29 electrode tip looks broken)
for dd=1:length(data.trial)
    tmpChan12 = data.trial{1,dd}(12,:);
    tmpChan29 = data.trial{1,dd}(29,:);
    data.trial{1,dd}(12,:) = tmpChan29;
    data.trial{1,dd}(29,:) = tmpChan12;
end


%%% resample the data (get 1 sequence long epoch)
% so we have an integer number of samples per cycle
% and define the trials (trl) based on the resampled data
cfg.trialdef.preStimDuration = 0; % required for resample_ssvep
cfg.newFs = 85*6; %Integer number of samples per monitor refresh (~500)
cfg.trialdef.epochLength = 1/experimentData(1).condInfo.stimTagFreq * 2; % size of the window for cutting trials (in seconds) = 2 cycles = 600ms
data = resample_only(cfg,data);

save([dataOut eegFiles(ff).name(1:end-4) '_needClean'],'data')

%%%%%%%%%%%%%%%%%%%
% artefact rejection
% first check for extrem channels/trials
cfg.layout = 'biosemi128.lay';
cfg.method = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg, data);

cfg.method = 'summary';
cfg.keepchannel = 'repair'; % had to modify ft_rejectvisual line 343 so that layout was taken into account
[data] = ft_rejectvisual(cfg, data); % if I want to change the way how the channels are interpolated then will have to do channel repair separately (will also not have to change the rejectvisual function)

% do the cleaning on 20 channels cap + 3 electrodes around the eyes
cfg.channel =  {'C28','C29','C30','C16','C17','C18','C4','C21','D4','D23','D19','A1','B22','B26','B4','A19','A7','A15','A23','A28','EXG2', 'EXG3','EXG4'}; % {'fp1','fp2','p4','fz','f3','t7','c3','cz','c4','t8','p4','pz','p3','o1','oz','o2'}
cfg.viewmode = 'vertical';
cfg = ft_databrowser(cfg,data);
cfg.artfctdef.reject = 'complete';
[cleanData] = ft_rejectartifact(cfg, data);

save([dataOut eegFiles(ff).name(1:end-4) '_clean'],'cleanData')

end



