function [data, trl] = resample_ssvep(cfg, datain)

trl = [];
preStimDuration = cfg.trialdef.preStimDuration;
newFs = cfg.newFs;
epochLengthSecs = cfg.trialdef.epochLength; % 2 s
trialLengthSecs = cfg.trialdef.trialLength; % 14 s
offset = 0;
cycleLength = cfg.trialdef.cycleLength*newFs;

% resample the entire trial
for iTrial = 1:length(datain.trial)
    newTrialLength = round(newFs*trialLengthSecs);
    newTimeBase = linspace(0,trialLengthSecs,newTrialLength+1);
    newTimeBase = newTimeBase(1:end-1);
    [time{iTrial}] = deal(newTimeBase);
end
resampleCfg.time = time;
resampleCfg.method = 'linear';
dataRes = ft_resampledata(resampleCfg, datain);

% redefine the trl (2 s epoch) based on the resampled trial 
% first remove the pre and post stimulus
numCyclesPerEpoch = (trialLengthSecs-preStimDuration*2) / epochLengthSecs;
for iTrial = 1:length(datain.trial)
    condition = datain.trialinfo(iTrial);
    for iCycle = 1:numCyclesPerEpoch
        begsample = round((preStimDuration*newFs+1)+epochLengthSecs*newFs*(iCycle-1)) + trialLengthSecs*newFs*(iTrial-1);
        endsample = begsample + epochLengthSecs*newFs-1;
        trl(end+1, :) = [begsample endsample offset condition iTrial cycleLength iCycle];
    end
end
        
cfg.trl = trl;
data = ft_redefinetrial(cfg,dataRes);