%%% create a trigger channel to test the SSVEP analysis
% you should not re-reference the data when doing so

%%%% pff so silly!! I already have it, it is channel 137!
%%% it's only during the SSVEP computation that I remove the status
%%% (trigger) channel. I could run it separately and see what SSVEP
%%% it gets

ff=18;
dataDir = '/Users/marleneponcet/Documents/data/MAE/originalData/';
dataOut = '/Users/marleneponcet/Documents/data/MAE/cleanData/';
eegFiles = dir([dataDir '*.bdf']);
behavFiles = dir([dataDir '*.mat']);

cfg.dataset   =  [dataDir eegFiles(ff).name];

hdr = ft_read_header(cfg.dataset);
schan = find(strcmpi(hdr.label,'STATUS'));

endsample = hdr.nSamples*hdr.nTrials;

sdata = ft_read_data(cfg.dataset, 'header', hdr, 'dataformat', 'biosemi_bdf', 'begsample', 1, 'endsample', endsample, 'chanindx', schan);

sdata = bitand(int32(sdata), 2^24-1);
byte1 = 2^8  - 1;
byte2 = 2^16 - 1 - byte1;
trigger = bitand(sdata, bitor(byte1, byte2)); % this is contained in the lower two bytes

trigValue = 1; % this is how it was set up for the experiment 
trigTest = bitand(sdata,trigValue);
trigTest = double(trigTest);



% read the behavioural file to get parameters for SSVEP
load([dataDir behavFiles(ff).name])
cfg.trialdef.freqTag = experimentData(1).condInfo.testFreq; % tagging freq (4.25 Hz = 85/20)
cfg.trialdef.cycleLength = 1/cfg.trialdef.freqTag; % duration one cycle (0.2353)
cfg.trialdef.trialLength = cfg.trialdef.cycleLength*21;
cfg.abortTrigger = 99; % invalid trial

% define trials
cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 112]; % all possible conditions
cfg.trialdef.ssvepTagVal = 1;
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'df_MAE';
[cfg] = ft_definetrial(cfg);
    
% pre-processing
cfg.reref         = 'no';
cfg.channel    = 1:10; 
cfg.lpfilter  = 'no';
cfg.lpfreq = 85; 
cfg.hpfreq = 1;
cfg.hpfilter = 'no'; 
cfg.detrend = 'no';
[data] = ft_preprocessing(cfg, trigTest);

cfg.channel = 137;
% cfg.channel = 'all';
[data] = ft_preprocessing(cfg);

