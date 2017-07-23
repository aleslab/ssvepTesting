clear all;
% To facilitate data-handling and distributed computing you can use
cfg.dataset   =  '/Users/ales/data/ssvepTesting/SSVEPtest3.bdf'

%cfg.dataset   =  '/Users/ales/data/ssvepTesting/testTrigger4/testTrigger4.bdf'
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.padding      = length (in seconds) to which the trials are padded for filtering (default = 0)
%   cfg.padtype      = string, type of padding (default: 'data' padding or
%                      'mirror', depending on feasibility)
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%




% cfg.demean        ='yes';
% cfg.reref         = 'yes';
% cfg.refchannel    = {'A1'};
cfg.lowpassfilter = 'no';
cfg.lpfreq        = 100;
cfg.demean        ='no';
cfg.reref         = 'no';

[data] = ft_preprocessing(cfg)

cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialfun = 'trialfun_ssvep'; 
[cfg] = ft_definetrial(cfg)

trl = cfg.trl;
cfg = [];
cfg.trl = trl;
data = ft_redefinetrial(cfg,data);

eegChan = ft_channelselection('A*',data.label);
nChan = length(data.label)
nEegChan = length(eegChan)

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

cfg.newFs = 85*6; %Integer number of samples per monitor refresh
%TODO: add code to read these values from the sesion info.
%This makes each cycle take an integer value;
data = resample_steadystate(cfg,data)
dataCycle = data;
[timelockCycle] = ft_timelockanalysis(cfg, data)

%These should be better designed. Ideally these should get picked when the
%paradigm is created. 
%For the contrast reversal exp: we had 57 cycles at 4.722 hz (.211)
%
cfg.trimdur = 6*.2117; %nycles to trim 
cfg.epochdur   = 9*.2117; %For an epoch duration of 1.9

data = create_epochs(cfg,data)

cfg.vartrllength = 2;
[timelockEpoch] = ft_timelockanalysis(cfg, data)

%remove the extra channels. 
cfg.channel =  {'all','-GSR1', '-GSR2', '-Erg1','-Erg2','-Resp','-Plet','-Temp', '-Status'};
dataEeg = ft_selectdata(cfg,data);
%Do the steadystate analysis
[Axx] = ft_steadystateanalysis(cfg, dataEeg)


%Now do a quick spectrum of the 12th channel:
figure(101);
pdSpecPlot(Axx.freq(2:80),Axx.Amp(12,2:80)',Axx.tCircPval(12,2:80)<.05)
title('Channel 12')

%Now setup an interactive plot
cfg.layout = 'biosemi32.lay';
cfg.channel = {'all'};
figure(201);clf;
interactiveTopoSpecPlot(cfg,Axx)


%Now lets look at the photodiode:
figure(102);
cfg.channel = { 'Erg1'}
dataPhotodiode = ft_selectdata(cfg,data);
[AxxPhotoDiode] = ft_steadystateanalysis(cfg, dataPhotodiode)
pdSpecPlot(AxxPhotoDiode.freq(2:80),AxxPhotoDiode.Amp(1,2:80)',AxxPhotoDiode.tCircPval(1,2:80)<.05)
title('PhotoDiode')

