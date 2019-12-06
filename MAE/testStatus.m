%%%% run analysis to test status channel



addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults


ff=2;

% path to the files
% dataDir = 'C:\Users\Marlene\Documents\dataStAndrews\MAE\';
dataDir = '/Users/marleneponcet/Documents/data/MAE/';
dataOut = '/Users/marleneponcet/Documents/data/MAE/cleanData/';
eegFiles = dir([dataDir '*.bdf']);
behavFiles = dir([dataDir '*.mat']);


cfg.dataset   =  [dataDir eegFiles(ff).name];
    
% read the behavioural file to get parameters for SSVEP
% load([dataDir behavFiles(ff).name])
load('/Users/marleneponcet/Documents/data/MAE/originalData/MAE_18_all_20191018.mat')
cfg.monitorRefresh = sessionInfo.expInfo.monRefresh;
cfg.framePerCycle = 20;% number of frames in one cycle
cfg.trialdef.freqTag = cfg.monitorRefresh/cfg.framePerCycle; % 4.25 Hz
cfg.trialdef.cycleLength = 1/cfg.trialdef.freqTag; % duration one cycle (0.2353)
cfg.trialdef.trialLength = cfg.trialdef.cycleLength*21; % 21 test cycles
% can potentially use the real duration of the trial but if so, also need
% to change the way the resampling is done by using the screen frequency
% and not a rounded one (e.g. not 85*6) 
cfg.abortTrigger = 99; % invalid trial

% Used for S15 to get the condition number
% tt = struct2table(experimentData);
% aa = struct2table(tt.condInfo);
% cfg.condition = aa.triggerCond';

% define trials
cfg.samplingRate = 2048;
cfg.trialdef.bitmask = 2^9-1; %Values to keep.
% cfg.trialdef.condRange = [101 112]; % all possible conditions
cfg.trialdef.condRange = [120 132]; % all possible conditions
cfg.trialdef.ssvepTagVal = 1;
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'df_MAE';
[cfg] = ft_definetrial(cfg);
    
% pre-processing
cfg.demean        ='no'; % useless with detrend. applies baseline correction to the entire trial by default (or to the specified window in cfg.baselinewindow)
cfg.reref         = 'no';
cfg.refchannel    = {'A1','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'}; % A3 = CPz / use Cz = A1
% cfg.refchannel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
cfg.lpfilter  = 'no';
cfg.lpfreq = 85; % screen 85Hz ...  49 would really clear noise, not 85
cfg.hpfreq = 1;
cfg.hpfilter = 'no'; % does NOT work anyway
cfg.detrend = 'no';
[data] = ft_preprocessing(cfg);

%%%%%%%%%%%%%% At the end of the preprocessing, the timeline for a trial is
%%%%%%%%%%%%%% up to 4.8240 s even though I defined it as 4.9398 s (see
%%%%%%%%%%%%%% trialLength value). Because of this, the last cycle cannot
%%%%%%%%%%%%%% be included in the analysis: I'll get NaN! There is no data after 4.8 s! 
   


% resample the data so we have an integer number of samples per cycle
% and define the trials (trl) based on the resampled data
cfg.newFs = cfg.monitorRefresh*6; %Integer number of samples per monitor refresh (~500)
cfg.trialdef.epochLength = 1/cfg.monitorRefresh*20*6; % size of the window for cutting trials (in seconds)
% cfg.trialdef.epochLength = 0.01176*20*9;
cfg.trialdef.preStimDuration = 1/cfg.monitorRefresh*20*2 ; 
% cfg.trialdef.preStimDuration = 0.01176*20*2;
data = resample_ssvep(cfg,data);


cfg.channel =  {'Status'};
allcond = unique(data.trialinfo(:,1));
    
for cond=1:length(allcond)
    cfg.trials = find(data.trialinfo(:,1) == allcond(cond));
    [statusChannel] = ft_steadystateanalysis(cfg, data);
end

cfg.channel =  {'Erg1'};
allcond = unique(data.trialinfo(:,1));
    
for cond=1:length(allcond)
    cfg.trials = find(data.trialinfo(:,1) == allcond(cond));
    [ergoChannel] = ft_steadystateanalysis(cfg, data);
end

figure;plot(data.trial{1,1}(139,:))
cc= [data.trial{:}];
dd= reshape(cc,144,[],96);
figure;plot(squeeze(dd(139,:,[1:3:end])))
figure; plot(ergoChannel.wave)
