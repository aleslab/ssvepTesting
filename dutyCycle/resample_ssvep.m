function [data, trl] = resample_ssvep(cfg, datain)

trl = [];
preStimDuration = cfg.trialdef.preStimDuration;
newFs = cfg.newFs;
epochLengthSecs = cfg.trialdef.epochLength; % 1.88 s
trialLengthSecs = cfg.trialdef.trialLength; % 11.67 s
offset = 0;

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
numCyclesPerEpoch = round((trialLengthSecs-preStimDuration*2) / epochLengthSecs); 
for iTrial = 1:length(datain.trial)
    condition = datain.trialinfo(iTrial,1);
    cycleLength = round(trialLengthSecs/datain.trialinfo(iTrial,4)*newFs); % round just so that it gives an integer number
    for iCycle = 1:numCyclesPerEpoch
        begsample = round((preStimDuration*newFs+1)+epochLengthSecs*newFs*(iCycle-1)) + trialLengthSecs*newFs*(iTrial-1);
        endsample = begsample + epochLengthSecs*newFs-1;
        trl(end+1, :) = [round(begsample) round(endsample) offset condition iTrial cycleLength iCycle]; % round just so that it gives an integer number
    end
end
        
cfg.trl = trl;
data = ft_redefinetrial(cfg,dataRes);