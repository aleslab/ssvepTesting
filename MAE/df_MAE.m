function trl = df_MAE(cfg)
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
startTest = cfg.trialdef.testTrigger;

%First find all condition starts
idx = 1;
for i=1:length(event)
    
    if isempty(event(i).value)
        continue;
    end
    
    % all tests (tagged 150) are changed into the condition number
    % the tag for the condition number is deleted 
    thisMaskedVal = bitand(event(i).value,bitmask);
    if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2)
        condition = thisMaskedVal;
    end
    if thisMaskedVal == startTest
        condNum(idx) = condition;
        condEventIdx(idx) = i;
        idx = idx+1;
    end
    
end

% there should be a total of 216 tests
condEventIdx(idx) = length(event); % create a last fake event (necessary to get the last trial included,see the following loop)

%Now lets break up each condition trial and find all ssvep cycles.
for iTrial = 1:length(condEventIdx)-1
    
    startCond = condEventIdx(iTrial);
    endCond = condEventIdx(iTrial+1); % iTrial+1: that is why there is a fake last event created earlier

    condValues  = [event(startCond:endCond).value];
    condValues  = bitand(condValues,bitmask); % trigger value
    condSamples = [event(startCond:endCond).sample]; % time of sampling
    
    if length(condSamples) ~= 22 && length(condSamples) ~= 25
        fprintf('test %d with wrong number of triggers \n ',iTrial)
    end        
    cycleStarts = condSamples(find(condValues==ssvepTagVal));
    cycleEnd = condSamples(21);
    
    if find(abs(diff(diff(cycleStarts))) > 5)   
        fprintf('test %d with wrong number of frames \n ',iTrial)
    end
    
    begsample = cycleStarts(1); % first trigger
    endsample = cycleEnd; % end of stim
    offset = 0;
    cycleLengthSamp = ceil(endsample - begsample);
    nbCycleStarts = length(cycleStarts);
    trl(end+1, :) = [begsample endsample offset condNum(iTrial) iTrial cycleLengthSamp nbCycleStarts];
    
end

