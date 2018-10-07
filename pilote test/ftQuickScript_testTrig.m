

clear all


addpath /home/marlene/Git/fieldtrip
ft_defaults
cfg.dataset = '/home/marlene/Documents/testTrigger4.bdf';


cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.trialfun = 'lock2SsvepTag_testTrig'; 

[cfg] = ft_definetrial(cfg);

cfg.demean        ='yes';
cfg.reref         = 'yes';
cfg.refchannel    = {'A32'}; % 32=Cz?
[data] = ft_preprocessing(cfg);

plot(data.trial{1,1}(36,:))
