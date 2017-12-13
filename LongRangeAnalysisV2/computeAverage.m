
% S11 is missing the "elec column"
% need to make sure that all Axx have the same fields

clear all;
addpath /Users/marleneponcet/Documents/Git/fieldtrip20170924
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/AxxFiles/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

for cond=1:6
    gp(cond) = doAverage(dataDir,listData,cond);
end
interactiveSteadyStatePlot(cfg,gp)

