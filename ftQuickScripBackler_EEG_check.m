clear cfg data;

% restoredefaultpath;                  % Clear the cluttered path

current_dir = fileparts(mfilename("fullpath"));

if isempty(current_dir)
    current_dir = pwd;
end



addpath("/Users/ales/git/fieldtrip-aleslab-fork");
addpath(genpath(current_dir));
ft_defaults;

% To facilitate data-handling and distributed computing you can use
% /Users/ales/git/ssvepTesting/data/test_alpha_with_displayed.bdf
% cfg.dataset   =  fullfile(current_dir,'data','test_audio_triggers.bdf')
% cfg.dataset   =  fullfile(current_dir,'data','test_alpha_with_displayed.bdf')

datadir = '/Users/ales/projects/backler_ssvep/03 EEG PILOT/01 Raw Data/EEG/'

cfg.dataset   =  fullfile(datadir,'av_sync_rse_S01_20261009.bdf');

%These paramaters are going to be applied to the raw continous data
%before clipping out trials. 

% data.trial{1} = abs(data.trial{1})




%Next were going to cut out the trials. Doing this in a second step to
%allow doing some processing before cutting out the trials. 
%First we're going to cut out all the individual cycles of the stimulus:
cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1; %event value to use to video ssvep cycles
% cfg.trialdef.ssvepTagVal = 2; %event value to use to audio tag ssvep cycles
% This is a link to a function that does the "trial" chopping and
% organizing for fieldtrip. 
cfg.trialfun = 'trialfun_ssvep'; 

%First define the trials:
[cfg] = ft_definetrial(cfg);
%Now take the trials out of the cfg structure
trl = cfg.trl; 


cfg.demean        ='no'; % useless with detrend. applies baseline correction to the entire trial by default (or to the specified window in cfg.baselinewindow)
cfg.reref         = 'yes';
cfg.refchannel    = {'A1','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'}; % A3 = CPz / use Cz = A1
cfg.lpfilter  = 'yes';
cfg.lpfreq = 85; % screen 85Hz ...  49 would really clear noise, not 85
cfg.hpfreq = 1;
cfg.hpfilter = 'no'; % does NOT work anyway
cfg.detrend = 'yes';

[data] = ft_preprocessing(cfg); %Load the data and do the above steps. 



%Need to reset the config for redefining the trials 
%I don't remember why. 
cfg = [];
cfg.trl = trl;
%And finally cut the trials out of the already loaded data. 
data = ft_redefinetrial(cfg,data);

%Now lets figure out how many EEG channels there are. We will do this by
%taking out the extra channels we saved. 
%This is kinda dumb.  We could grab jjust the "EEG" channels,
%But For the trigger tests I want to analyze the Erg input as well as the
%EEG
chanSel =  {'eeg'};%{'all','-GSR1', '-GSR2','-Erg2','-Resp','-Plet','-Temp', '-Status'};
eegChan = ft_channelselection(chanSel,data.label);
nChan = length(data.label);
nEegChan = length(eegChan);


%In order to keep the extra channels in we need to make a montage that
%doesn't apply the reference to the extra channels. 
%This is only needed for EEG processing so we'll comment it out for the
%trigger stuff
%aveRefMtx   = eye(nEegChan)-ones(nEegChan)/nEegChan;
% czRefMtx      = eye(nEegChan);
% czRefMtx(:,32) =czRefMtx(:,32)-1; 
% otherMtx    = eye(nChan-nEegChan);
% montage.tra = blkdiag(czRefMtx,otherMtx);
% montage.labelold = data.label;
% montage.labelnew = data.label;
% data = ft_apply_montage(data,montage);

%Next resample the data so we have an integer number of samples per cycle.

%Lets read it from the data file 
ptbCorgiInfo = load("/Users/ales/projects/backler_ssvep/03 EEG PILOT/01 Raw Data/Behaviour/P1/avEEGpilot_S001_expPilot_20261009_115133.mat")

% ptbCorgiInfo = load("/Users/ales/projects/backler_ssvep/03 EEG PILOT/01 Raw Data/Behaviour/P2/avEEGpilot_S002_expPilot_20261009_150832.mat");

actualFrameRate = ptbCorgiInfo.sessionInfo.expInfo.monRefresh;


cfg.newFs = actualFrameRate*6; %Integer number of samples per monitor refresh
%TODO: add code to read these values from the sesion info.
%We are going to resample the data to be an integer number of samples per
%cycle with a fixed time from the cycle start tag. 
data = resample_steadystate(cfg,data);


%Now that we have cycles with an equal number of samples we can combine the
%cycles together to form "epochs", an epoch forms the base time window used
%for the fourier transform.  These are 

%The values beow should be better designed. Ideally these should get picked when the
%paradigm is created. Gotta be careful here that everything works out 

cfg.trimdur = 5*(24/actualFrameRate); 
cfg.epochdur   = 10*(24/actualFrameRate); 

% We can calculate 
% 
data = create_epochs(cfg,data);


%DO EPOCH ARTIFACT REJECTION HERE

% cfg.vartrllength = 2;
% [timelockEpoch] = ft_timelockanalysis(cfg, data);

%remove the extra channels. 
cfg.channel =  {'eeg'};%{'all','-unkown','-GSR1', '-GSR2','-Erg2','-Resp','-Plet','-Temp', '-Status'};
dataEeg = ft_selectdata(cfg,data);
%Do the steadystate analysis
cfg.f1hz = actualFrameRate/24; %TODO: pull this from session files. 
cfg.f2hz = actualFrameRate/30;
% cfg.f1hz = 75/30; %TODO: pull this from session files. 

allcond = unique(dataEeg.trialinfo(:,1));
for cond=1:length(allcond)
    cfg.trials = find(dataEeg.trialinfo(:,1) == allcond(cond));
    [Axx(cond)] = ft_steadystateanalysis(cfg, dataEeg); % the field elec is present only if a channel has been replaced
end




%Now setup an interactive plot
%The trigger testing datasets didn't actually have electrodes plugged in. 
%We did plug the microphone into the erg1 channel for testing
%SO we're screwing up the layout info and spoofing stuff here just to get
%the plotting function to work.  Only the last channel has any real data. 
% cfg.layout = 'biosemi64.lay';
% cfg.channel = 1:Axx.nchan;
% figure;clf;
% 
% % There's a stupid revision here with a dumb 2 name.   We'll fix this soon
% interactiveSteadyStatePlot2(cfg,Axx)
% 
% %Now lets look at the photodiode:
% figure;
% cfgErgo= cfg;
% cfgErgo.channel = { 'Erg1'}
% dataPhotodiode = ft_selectdata(cfgErgo,data);
% [AxxPhotoDiode] = ft_steadystateanalysis(cfgErgo, dataPhotodiode)
% subplot(3,1,1);
% pdSpecPlot(AxxPhotoDiode.freq(2:80),AxxPhotoDiode.amp(1,2:80)',AxxPhotoDiode.pval(1,2:80)<.05)
% title('PhotoDiode')
% subplot(3,1,2);
% pdSpecPlot(AxxPhotoDiode.freq(80:200),AxxPhotoDiode.amp(1,80:200)',AxxPhotoDiode.pval(1,80:200)<.05);
% subplot(3,1,3);
% plot(AxxPhotoDiode.time,AxxPhotoDiode.wave)

