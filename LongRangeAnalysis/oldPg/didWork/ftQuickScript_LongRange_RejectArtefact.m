clear cfg data;
% cfg.dataset   =  'D:\StAndrews\LongRangeSSVEP\LongRangeS01.bdf';
% addpath D:\GitHubRepo\fieldtrip
% addpath D:\GitHubRepo\ssvepTesting\svndlCopy
% addpath D:\GitHubRepo\ssvepTesting\biosemiUpdated

% noisy trials: 280, 284, 394, 506, 586, 616 

addpath /Users/marleneponcet/Documents/Git/fieldtrip
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
cfg.dataset   =  '/Users/marleneponcet/Documents/LongRangeSSVEP/LongRangeS01.bdf';

ft_defaults

% read the behavioural file
load('longRange_01__20170727_153334.mat')
tableData = struct2table(experimentData);

% save some of the data from the behavioural file into the cfg
cfg.trialdef.preStimDuration = experimentData(1).trialData.preStimDuration;
cfg.trialdef.nbTotalCycles = experimentData(1).trialData.nbTotalCycles;

% define trials 
cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'lock2ssvep_LongRange_V2'; 
[cfg] = ft_definetrial(cfg);

% pre-processing
cfg.lpfilter  = 'yes';
cfg.lpfreq        = 80;
cfg.demean        ='yes';
cfg.channel    = 1:128;
cfg.reref         = 'yes'; 
cfg.refchannel    = {'A1'}; % A3 = CPz / use Cz = A1
[data] = ft_preprocessing(cfg); 

% resample the data so we have an integer number of samples per cycle.
cfg.newFs = 85*8; %Integer number of samples per monitor refresh
data = resample_steadystate_test(cfg,data);
% because of reseampling, should change data.trialinfo
data.trialinfo(:,3) = size(data.time{1,1},2) / 5; % 5 cycles / epoch


%DO EPOCH ARTIFACT REJECTION 
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
cfg.artfctdef.eog.channel = 'C17';
cfg.artfctdef.eog.feedback = 'yes';
% cfg.artfctdef.eog.bpfreq     = [1 15];
% cfg.artfctdef.eog.bpfiltord  = 4;
[cfg, artifact] = ft_artifact_eog(cfg,data);

cfg.artfctdef.reject = 'complete';
[cleanData] = ft_rejectartifact(cfg, data);

% optional: last look
cfg = ft_databrowser(cfg,cleanData); 
cfg.artfctdef.reject  = 'complete';
cleanData = ft_rejectartifact(cfg,cleanData); 

save('cleanData_S01', 'cleanData')


% cfg = [];
% cfg.vartrllength = 2;
% [timelockEpoch] = ft_timelockanalysis(cfg, cleandata);


% do the steadystate analysis for each conditions
% cfg.trimdur = 0;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
allcond = unique(cleanData.trialinfo(:,1));
for cond=1:length(allcond)
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
    [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData);
end
save('AxxS01', 'Axx')

cfg = [];
cfg.layout = 'biosemi128.lay';
interactiveSteadyStatePlot(cfg,Axx)
