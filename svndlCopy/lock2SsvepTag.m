function trl = lock2SsvepTag(cfg);
 
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
 
trl = [];
 
%Silly hard coded numbers
%
bitmask = 2^9-1;%Values to keep.
condRange = [101 165];
ssvepTagVal = 1;
epochLength = 2; %approx epoch length in seconds. 



%First find all condition starts
idx = 1;
for i=1:length(event)

    if isempty(event(i).value);
        continue;
    end
    
    thisMaskedVal = bitand(event(i).value,bitmask);
    
    if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2);
        condEventIdx(idx) = i;
        condNum(idx)      = thisMaskedVal;
        idx = idx+1;
    end
    
    
end

condEventIdx(idx) = length(event)+1;

%Now lets break up each condition trial and find all ssvep cycles. 
for iTrial = 1:length(condEventIdx)-1,

    startCond = condEventIdx(iTrial)+1;
    endCond = condEventIdx(iTrial+1)-1;
    
     condValues  = [event(startCond:endCond).value];
     condValues  = bitand(condValues,bitmask);
     condSamples = [event(startCond:endCond).sample];
    
    trial(iTrial).cycleStarts = condSamples(condValues==ssvepTagVal);
    trial(iTrial).condNum = condNum(iTrial);
    
    numCycles = length(trial(iTrial).cycleStarts);
    intMultiples = divisors(numCycles);
    %Length of time SSVEP is on:
    ssvepTime = (trial(iTrial).cycleStarts(end)-trial(iTrial).cycleStarts(1))/hdr.Fs;
    numEpochsByTime = ssvepTime/epochLength;
    
    intIdx = find((intMultiples-numEpochsByTime)>=0,1,'first');
    numCyclesPerEpoch = intMultiples(intIdx)
    
    meanSampPerCycle= round(mean(diff(trial(iTrial).cycleStarts)));
    
    for iCycle = 1:numCyclesPerEpoch:numCycles;
        
        
        begsample     = trial(iTrial).cycleStarts(iCycle);
        endsample     = begsample+numCyclesPerEpoch*meanSampPerCycle;
        offset        = 0;
        trl(end+1, :) = round([begsample endsample offset]);
    end
    
    %Now find how many epochs fit into the data:
    
end

