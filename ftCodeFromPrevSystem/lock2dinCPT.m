function [trl, eventName] = lock2dinCPT(cfg);
 
hdr   = read_header(cfg.dataset);
event = read_event(cfg.dataset);
 
trl = [];
eventName = [];
 
begTime = cfg.trialdef.prestim;
stopTime = cfg.trialdef.poststim;
%numConds = length(cfg.trialdef.eventvalue);
%Conds = cfg.trialdef.eventvalue;

for i=1:length(event)
    if strcmp(event(i).type, 'trigger')
        % it is a trigger, see whether it has the right value
        if any(strcmp(event(i).value, cfg.trialdef.eventvalue))

            if strcmp(event(i+1).value,'DIN4');

                % add this to the trl definition


				eventName{i}  = event(i).value;
                begsample     = event(i+1).sample - begTime*hdr.Fs;
                endsample     = event(i+1).sample + stopTime*hdr.Fs - 1;
                offset        = -begTime*hdr.Fs;
                trl(end+1, :) = round([begsample endsample offset]);
            else
                error('DIN4 does not follow condition trigger');
            end
        end

    end
end

results = regexp(eventName, 'c\d\d\d');
resIdx = ~cellfun('isempty', results);


eventName = eventName(resIdx);
