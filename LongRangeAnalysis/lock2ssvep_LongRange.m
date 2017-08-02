function trl = lock2ssvep_LongRange(cfg)

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

trl = [];

%Silly hard coded numbers
%
bitmask = 2^9-1;%Values to keep.
condRange = [101 165];
ssvepTagVal = cfg.trialdef.ssvepTagVal;
epochLength = cfg.trialdef.epochLength;
preStimDuration = cfg.trialdef.preStimDuration;
nbTotalCycles = cfg.trialdef.nbTotalCycles;

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
    startCond = condEventIdx(iTrial)+1;
    endCond = condEventIdx(iTrial+1)-1; % iTrial+1: that is why there is a fake last event created earlier
    
    condValues  = [event(startCond:endCond).value];
    condValues  = bitand(condValues,bitmask);
    
    if find(condValues == 98) % check if it is any invalid trial
        fprintf('trial %d invalid \n',iTrial)
    else
        indexTrial = find(condValues == 64);
        condSamples = [event(startCond:endCond).sample];
        cycleStarts = condSamples(find(condValues==ssvepTagVal));
        
        if length(indexTrial) ~= 2 || condValues(indexTrial(2)-1) ~= 10
            fprintf('\n trial %d with missing start or end trigger',iTrial)
        elseif find(abs(diff(diff(cycleStarts))) > 5) % check for a missing frame
            fprintf('\n trial %d with a missing frame',iTrial)
        elseif length (cycleStarts) ~= nbTotalCycles
            fprintf('\n not the expected number of cycles in trial %d',iTrial)
        else % no problem with this trial
            meanSampPerCycle= round(mean(diff(cycleStarts)));
            % only extract epochs within the pre and post stimulus time
            nbSkipCycles = round (preStimDuration / (meanSampPerCycle/hdr.Fs));
            stimCycles = cycleStarts(nbSkipCycles+1 : nbTotalCycles - nbSkipCycles);
%             stimCyclesPlusOne = cycleStarts(nbSkipCycles+1 : nbTotalCycles - nbSkipCycles + 1); % add the start of next sample so that it can be used in the next iCycle loop
            
            numCyclesPerEpoch = floor(epochLength / (meanSampPerCycle/hdr.Fs));
            
            trial(numTrial).cycleStarts = stimCycles;
            trial(numTrial).condNum = condNum(iTrial);
            
            for iCycle = 1:numCyclesPerEpoch:length(stimCycles)
                begsample     = trial(numTrial).cycleStarts(iCycle);
                endsample = begsample+meanSampPerCycle*numCyclesPerEpoch-1;
%                 endsample     = stimCyclesPlusOne(iCycle+1)-1; %One sample before the next cycle.
                offset        = 0;
                cycleLengthSamp = endsample - begsample +1; % have to add 1 because the cycle is with +1 (since remove the one above? - still confused)
                trl(end+1, :) = [begsample endsample offset trial(numTrial).condNum iTrial cycleLengthSamp iCycle];
            end
            numTrial = numTrial + 1;
        end
    end
end

