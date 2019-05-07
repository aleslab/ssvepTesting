function [data, trl] = resample_only(cfg, datain)
% only resample data without redefining trial

trl = [];
preStimDuration = cfg.trialdef.preStimDuration;
newFs = cfg.newFs;
trialLengthSecs = cfg.trialdef.epochLength; 
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
data = ft_resampledata(resampleCfg, datain);
