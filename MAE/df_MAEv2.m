function trl = df_MAEv2(cfg)
% v2: use the last SSVEP trigger as a "end of cycle" trigger

% condition number is at the beginning of the trial (101-112)
% 150 is a trigger sent at the beggining of each test 
% so for one trial there should be a condition number followed by 9x150 
% (there are 9 tests)

event = ft_read_event(cfg.dataset);
trl = [];

%Silly hard coded numbers
bitmask = 2^9-1;%Values to keep.
condRange = cfg.trialdef.condRange;
ssvepTagVal = cfg.trialdef.ssvepTagVal;
abortExpTrigger = cfg.abortTrigger;
frameDuration = 1/cfg.monitorRefresh;
samplingRate = 1/cfg.samplingRate;
samplesPerFrame = round(frameDuration/samplingRate);


% use frameDuration to compute nb of samples in one frame

%First find all condition starts
idx = 1;
for i=1:length(event)
    
    %%%%% the trigger from Bits is read only at the end of the frame, that
    %%%%% is, it is delayed by one frame compared to the stimulus
    %%%%% presentation. To get it in sync, remove the number of samples
    %%%%% corresponding to one frame
    event(i).sample = event(i).sample - samplesPerFrame;  
    
    if isempty(event(i).value)
        continue;
    end
    
    thisMaskedVal = bitand(event(i).value,bitmask);
%     if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2)
%         condNum(idx) = thisMaskedVal;
%         condEventIdx(idx) = i;
%         idx = idx+1;
%     end
    if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2)
        prevCond = bitand(event(i-1).value,bitmask);
%         condNum(idx) = thisMaskedVal + prevCond-100;
        condNum(idx) = thisMaskedVal;
        condEventIdx(idx) = i;
        idx = idx+1;
    end
    
end

% condNum = condNum+cfg.condition;


% there should be a total of 216 tests
condEventIdx(idx) = length(event) +1; % create a last fake event (necessary to get the last trial included,see the following loop)

%Now lets break up each condition trial and find all ssvep cycles.
for iTrial = 1:length(condEventIdx)-1
    
    startCond = condEventIdx(iTrial);
    endCond = condEventIdx(iTrial+1)-1; % iTrial+1: that is why there is a fake last event created earlier
    % and -1 because event(created trial).value doesn't exist

    condValues  = [event(startCond:endCond).value];
    condValues  = bitand(condValues,bitmask); % trigger value

    if find(condValues == abortExpTrigger) % check if it is any invalid trial
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
        
        if length(indexTrial) ~= 2 
            fprintf('trial %d with missing start or end trigger \n',iTrial)
        elseif find(abs(diff(diff(cycleStarts))) > 5) % check for a missing frame
            fprintf('trial %d with a missing frame \n ',iTrial)
        else % no problem with this trial
            begsample = condSamples(indexTrial(1) + 1); % first trigger
            endsample = condSamples(indexTrial(2)-1); % end of stim
%             endsample = condSamples(indexTrial(2)-1); % end of stim 
            offset = 0;
            cycleLengthSamp = ceil(endsample - begsample);
            nbCycleStarts = length(cycleStarts)-1;
            trl(end+1, :) = [begsample endsample offset condNum(iTrial) iTrial cycleLengthSamp nbCycleStarts];
        end        
        
    end
    
end

