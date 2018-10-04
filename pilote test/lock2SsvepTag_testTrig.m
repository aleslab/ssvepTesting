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
    allVal(i) = thisMaskedVal;
    
    if thisMaskedVal>=condRange(1) && thisMaskedVal<=condRange(2)
        condEventIdx(idx) = i;
        condNum(idx)      = thisMaskedVal;
        idx = idx+1;
    end
    
end
condEventIdx(idx) = length(event)+1; % create a last fake event (necessary to get the last trial included,see the following loop)
    
    condValues  = [event(1:end).value]; 
    condValues  = bitand(condValues,bitmask);
    condSamples = [event(1:end).sample];
    diff(condSamples)/hdr.Fs

trl(1, :) = round([condSamples(1) condSamples(end) 0]);

end

