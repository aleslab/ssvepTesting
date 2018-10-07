function trl = lock2SsvepTag_LR(cfg)

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

trl = [];

%Silly hard coded numbers
%
bitmask = 2^9-1;%Values to keep.
condRange = [101 165];
ssvepTagVal = 1;
epochLength = 2; %approx epoch length in ms.
numExpectedCycles = 15;

% check if there is any invalid trial
for i=1:length(event)
    thisMaskedVal = bitand(event(i).value,bitmask);
    if thisMaskedVal == 98
        disp('invalid trials, do something')
    end
end

%First find all condition starts
idx = 1;
for i=1:length(event)
    
    if isempty(event(i).value)
        continue;
    end
    
    thisMaskedVal = bitand(event(i).value,bitmask);
    
    if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2)
        condEventIdx(idx) = i;
        condNum(idx)      = thisMaskedVal;
        idx = idx+1;
    end
    
end

condEventIdx(idx) = length(event)+1; % create a last fake event (necessary to get the last trial included,see the following loop)

numTrial = 1;
%Now lets break up each condition trial and find all ssvep cycles.
for iTrial = 1:length(condEventIdx)-1
    % get the ssvep triggers for this trial
    % +1/-1 should remove the condition number, not the experiment start/stop trigger.
    % Also -1 necessary for the last trial: the last event does not have a value (that's the fake event created just before the loop)
    startCond = condEventIdx(iTrial)+1;
    endCond = condEventIdx(iTrial+1)-1;
    
    condValues  = [event(startCond:endCond).value];
    condValues  = bitand(condValues,bitmask);
    condSamples = [event(startCond:endCond).sample];
    
    %     % get only the 'real' ssvep triggers
    %     cycleStarts = condSamples(find(condValues<5));  % include dim flashes
    %     numCycles = length(cycleStarts);
    %     %     trial(iTrial).cycleStarts = condSamples(condValues==ssvepTagVal); % does not include dim flashes
    %
    %     % sanity checks
    %     if  numCycles ~= numExpectedCycles % check for a missing cycle
    %         fprintf('\n not the expected number of cycles in trial %d',iTrial)
    %     elseif find(abs(diff(diff(cycleStarts))) > 5) % check for a missing frame
    %         fprintf('\n wrong duration trial %d',iTrial)
    %     else
    %         trial(numTrial).cycleStarts = cycleStarts;
    %         trial(numTrial).condNum = condNum(iTrial);
    %
    %         %     numCycles = length(trial(iTrial).cycleStarts);
    %         meanSampPerCycle= mean(diff(trial(numTrial).cycleStarts));
    %         numCyclesPerEpoch = floor(epochLength / (meanSampPerCycle/hdr.Fs));
    %
    %         for iCycle = 1:numCyclesPerEpoch:numCycles-(numCyclesPerEpoch-1)
    %             begsample     = trial(numTrial).cycleStarts(iCycle);
    %             endsample     = begsample+numCyclesPerEpoch*meanSampPerCycle;
    %             offset        = 0;
    %             trl(end+1, :) = round([begsample endsample offset trial(numTrial).condNum]);
    %         end
    %         numTrial = numTrial + 1;
    %     end
    
    % SHOULDN'T BE DONE THIS WAY
    cycleStarts = condSamples(find(condValues<5));
    firstTrigger = cycleStarts(1);
    trial(iTrial).condNum = condNum(iTrial);
    
    meanSampPerCycle= 819; % mean(diff(trial(numTrial).cycleStarts));
    numCyclesPerEpoch = floor(epochLength / (meanSampPerCycle/hdr.Fs));
    numCycles = 15;
    
    for iCycle = 1:numCyclesPerEpoch:numCycles-(numCyclesPerEpoch-1)
        begsample     = firstTrigger+(iCycle-1)*meanSampPerCycle;
        endsample     = begsample+numCyclesPerEpoch*meanSampPerCycle-1;
        offset        = 0;
        trl(end+1, :) = round([begsample endsample offset trial(iTrial).condNum iTrial meanSampPerCycle]);
    end
    
end

