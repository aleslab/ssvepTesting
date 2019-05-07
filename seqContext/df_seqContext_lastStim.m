function trl = df_seqContext_lastStim(cfg)

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
sbjNb = str2num(cfg.dataset(end-14:end-13));
trl = [];

bitmask = cfg.trialdef.bitmask;
condRange = cfg.trialdef.condRange;
ssvepTagVal = cfg.trialdef.ssvepTagVal;

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

%Now lets break up each condition trial and find all ssvep cycles.
for iTrial = 1:length(condEventIdx)-1
    % get the ssvep triggers for this trial +1/-1 should remove the
    % condition number, not the experiment start/stop trigger.
    startCond = condEventIdx(iTrial)+1;
    endCond = condEventIdx(iTrial+1)-1; % iTrial+1: that is why there is a fake last event created earlier

    condValues  = [event(startCond:endCond).value];
    condValues  = bitand(condValues,bitmask);
    
    if find(condValues == 98) % check if it is any invalid trial
        fprintf('trial %d invalid \n',iTrial)
    else
        indexTrial = find(condValues == 64);
        
        % should only include STATUS (no 'CM out of range') !!! 
        % condSamples = [event(startCond:endCond).sample];
        condSamples = [];
        for tt=startCond:endCond
            if strcmp(event(tt).type,'STATUS')
                condSamples = [condSamples, event(tt).sample];
            elseif isempty(condSamples)
                fprintf('trial %d CM out of range before trial \n ',iTrial)
                event(tt).type  
            else
                fprintf('trial %d CM out of range DURING trial \n ',iTrial)
                event(tt).type  
            end
        end
        cycleStarts = condSamples(find(condValues==ssvepTagVal));
        
        if length(indexTrial) ~= 2 || condValues(indexTrial(2)-1) ~= 10
            fprintf('trial %d with missing start or end trigger \n',iTrial)
        elseif find(abs(diff(diff(cycleStarts))) > 5) % check for a missing frame
            fprintf('trial %d with a missing frame \n ',iTrial)
        else % no problem with this trial
            % look for seq with 4 0
            stimLoc = [4 0];
            indS = strfind(cfg.stimSeq(iTrial,:),stimLoc);
            for locS=1:length(indS)
                begsample = condSamples(indexTrial(1) + indS(locS)); % indexTrial(1) is the trigger for trial start
                endsample = condSamples(indexTrial(1) + indS(locS) + 2); % strop 2 triggers later
                offset = 0;
                cycleLengthSamp = ceil(endsample - begsample);
                trl(end+1, :) = [begsample endsample offset condNum(iTrial) iTrial cycleLengthSamp locS];
            end
        end
    end
end

