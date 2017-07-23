function [trl event] = lock2SsvepTag(cfg);
 
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
 
trl = [];
 
%Silly hard coded numbers
%TODO: CHange these ot cfg.
bitmask      = ft_getopt(cfg.trialdef, 'bitmask',  2^9-1);
condRange    = ft_getopt(cfg.trialdef, 'condRange',   [101 165]);
ssvepTagVal = ft_getopt(cfg.trialdef, 'ssvepTagVal',   1);


%First mask the values. and replace them. Doing this in separate loop
%instead of a vectorized operation becaust event.value may be empty (e.g.
%CM_IN_RANGE events
for iEvent=1:length(event)
    if isempty(event(iEvent).value);
        continue;
    end
    event(iEvent).value = bitand(event(iEvent).value,bitmask);
end

%Now find how many events qualify as "conditions"
condEventIdx = [event.value]>= condRange(1) & [event.value]<= condRange(2); %Index to all events in cond range
condEventIdx = find(condEventIdx); %Convert form logical indexing to numeric indexing
condNumList = unique( [event(condEventIdx).value]);  %find unique conditions

%Now count how many "trials" exist for each condition
for iCond = 1:length(condNumList)
    trialsPerCond(iCond) = sum( [event(condEventIdx).value] == condNumList(iCond));
end

% %Now got through every condition
% idx = 1;
% for iEvent=1:length(event)
% 
%     if isempty(event(iEvent).value);
%         continue;
%     end
%     
%     
%     if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2);
%         condEventIdx(idx) = iEvent;
%         condNum(idx)      = thisMaskedVal;
%         idx = idx+1;
%     end
%     
%     
% end

%Silly adding 1 just because we don't have a condition end tag we'll use the end
%of the data as an implicit condition end.  
condEventIdx(end+1) = length(event)+1;

condTrialCount = zeros(length(condNumList),1); %Vector to hold the trial counts
%Now lets break up each condition trial and find all ssvep cycles. 
for iTrial = 1:length(condEventIdx)-1,

    startCond = condEventIdx(iTrial)+1;
    endCond = condEventIdx(iTrial+1)-1;
    
     condValues  = [event(startCond:endCond).value];    
     condSamples = [event(startCond:endCond).sample];
    
    trial(iTrial).cycleStarts = condSamples(condValues==ssvepTagVal);
    trial(iTrial).condNum = event(condEventIdx(iTrial)).value;
    condIdx = find(trial(iTrial).condNum ==condNumList); %Find which index this condition matches
    condTrialCount(condIdx) = condTrialCount(condIdx)+1; %increment trial counter by one
    
    numCycles = length(trial(iTrial).cycleStarts);
%     intMultiples = divisors(numCycles);
%     %Length of time SSVEP is on:
%     ssvepTime = (trial(iTrial).cycleStarts(end)-trial(iTrial).cycleStarts(1))/hdr.Fs;
%     numEpochsByTime = ssvepTime/epochLength;
%     
%     intIdx = find((intMultiples-numEpochsByTime)>=0,1,'first');
%     numCyclesPerEpoch = intMultiples(intIdx)
    
    meanSampPerCycle= round(mean(diff(trial(iTrial).cycleStarts)));
    
%    for iCycle = 1:numCyclesPerEpoch:numCycles;
    for iCycle = 1:numCycles;
        
        
        begsample     = trial(iTrial).cycleStarts(iCycle);
        
        if iCycle ==numCycles %On the last cycle start we don't have an end cycle so just use the mean length. 
            endsample     = begsample+meanSampPerCycle-1;%One sample before the next cycle.
        else
            endsample     = trial(iTrial).cycleStarts(iCycle+1)-1;%One sample before the next cycle.
        end
        
        offset        = 0;
        condNum       = trial(iTrial).condNum;
        trialNum      = condTrialCount(condIdx); 
        trl(end+1, :) = [begsample endsample offset condNum trialNum iCycle ];
    end
    
    
    %Now find how many epochs fit into the data:
    
end

