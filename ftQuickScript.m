clear cfg data;
% To facilitate data-handling and distributed computing you can use
cfg.dataset   =  '/Users/ales/data/ssvepTesting/SSVEPtest3.bdf'


%These paramaters are going to be applied to the raw continous data
%before clipping out trials. 
cfg.lowpassfilter = 'yes';
cfg.lpfreq        = 100;
cfg.demean        ='yes';
cfg.reref         = 'no'; 

[data] = ft_preprocessing(cfg); %Load the data and do the above steps. 


%Next were going to cut out the trials. Doing this in a second step to
%allow doing some processing before cutting out the trials. 
%First we're going to cut out all the individual cycles of the stimulus:
cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1; %event value to use to tag ssvep cycles
cfg.trialfun = 'trialfun_ssvep'; 

%First define the trials:
[cfg] = ft_definetrial(cfg);
%Now take the trials out of the cfg structure
trl = cfg.trl; 
cfg = [];
cfg.trl = trl;
%And finally cut the trials out of the already loaded data. 
data = ft_redefinetrial(cfg,data);

%Now lets figure out how many EEG channels there are. We will do this by
%taking out the extra channels we saved. 
chanSel =  {'all','-GSR1', '-GSR2', '-Erg1','-Erg2','-Resp','-Plet','-Temp', '-Status'};
eegChan = ft_channelselection(chanSel,data.label);
nChan = length(data.label);
nEegChan = length(eegChan);

%In order to keep the extra channels in we need to make a montage that
%doesn't apply the reference to the extra channels. 
%aveRefMtx   = eye(nEegChan)-ones(nEegChan)/nEegChan;
czRefMtx      = eye(nEegChan);
czRefMtx(:,32) =czRefMtx(:,32)-1; 
otherMtx    = eye(nChan-nEegChan);
montage.tra = blkdiag(czRefMtx,otherMtx);
montage.labelold = data.label;
montage.labelnew = data.label;
data = ft_apply_montage(data,montage);

%Next resample the data so we have an integer number of samples per cycle.

cfg.newFs = 85*8; %Integer number of samples per monitor refresh
%TODO: add code to read these values from the sesion info.
%We are going to resample the data to be an integer number of samples per
%cycle with a fixed time from the cycle start tag. 
data = resample_steadystate(cfg,data);


%Now that we have cycles with an equal number of samples we can combine the
%cycles together to form "epochs", an epoch forms the base time window used
%for the fourier transform.  These are 

%The values beow should be better designed. Ideally these should get picked when the
%paradigm is created. 
%For the contrast reversal exp: we had 57 cycles at 4.722 hz (.211). That
%means if we trim 6 cycles from the front and back we're left with 45
%cycles that divide neatly into 9 cycle epochs. 
%
cfg.trimdur = 6*.2117; %Trim 6 cycles from front and back.
cfg.epochdur   = 9*.2117; %For an epoch duration of 1.9

data = create_epochs(cfg,data)


%DO EPOCH ARTIFACT REJECTION HERE

cfg.vartrllength = 2;
[timelockEpoch] = ft_timelockanalysis(cfg, data)

%remove the extra channels. 
cfg.channel =  {'all','-GSR1', '-GSR2', '-Erg1','-Erg2','-Resp','-Plet','-Temp', '-Status'};
dataEeg = ft_selectdata(cfg,data);
%Do the steadystate analysis
[Axx] = ft_steadystateanalysis(cfg, dataEeg)


%Now do a quick spectrum of the 12th channel:
figure;
pdSpecPlot(Axx.freq(2:80),Axx.Amp(12,2:80)',Axx.tCircPval(12,2:80)<.05)
title('Channel 12')

%Now setup an interactive plot
cfg.layout = 'biosemi32.lay';
cfg.channel = {'all'};
figure;clf;
interactiveTopoSpecPlot(cfg,Axx)

interactiveSteadyStatePlot(cfg,Axx)

%Now lets look at the photodiode:
figure;
cfgErgo= cfg;
cfgErgo.channel = { 'Erg1'}
dataPhotodiode = ft_selectdata(cfgErgo,data);
[AxxPhotoDiode] = ft_steadystateanalysis(cfgErgo, dataPhotodiode)
subplot(3,1,1);
pdSpecPlot(AxxPhotoDiode.freq(2:80),AxxPhotoDiode.Amp(1,2:80)',AxxPhotoDiode.tCircPval(1,2:80)<.05)
title('PhotoDiode')
subplot(3,1,2);
pdSpecPlot(AxxPhotoDiode.freq(80:200),AxxPhotoDiode.Amp(1,80:200)',AxxPhotoDiode.tCircPval(1,80:200)<.05);
subplot(3,1,3);
plot(AxxPhotoDiode.time,AxxPhotoDiode.Wave)

