function [trl, eventName, reactionTime] = responseLockCPT(cfg);
 
hdr   = read_header(cfg.dataset);
event = read_event(cfg.dataset);
 
trl = [];
eventName = [];
timingCorrection = 432/1000;
 
%begTime = cfg.trialdef.prestim;
%stopTime = cfg.trialdef.poststim;
begTime = 1.5 * timingCorrection;
stopTime = .1;
%numConds = length(cfg.trialdef.eventvalue);
%Conds = cfg.trialdef.eventvalue;
a = 1;
for i=1:length(event)
    if strcmp(event(i).type, 'trigger')
        % it is a trigger, see whether it has the right value
        if any(strcmp(event(i).value, cfg.trialdef.eventvalue))

            if strcmp(event(i+2).value,'DIN1');

                % add this to the trl definition

                if strcmp(event(i+1).value,'DIN4');
                    
                    startTrial = event(i+1).sample;
                    buttonLatencySamples = event(i+2).sample - startTrial;
					%this is a bad way to do it (432 should be automated).
					%for CPT 432 is the number of samples in 1 second.
					%subtracting it here measures the button press in
					%milliseconds after target instead of milliseconds after
					%trial start
                    buttonLatencySeconds = (buttonLatencySamples - 432)/timingCorrection;
                    reactionTime(a).buttonLatencySamples = buttonLatencySamples;
                    reactionTime(a).buttonLatencySeconds = buttonLatencySeconds;
					a = a+1;
					
                else
                    disp('found a DIN1 with no preceeding DIN4');
				end
                
				tmpButtonLatencySamples = [reactionTime.buttonLatencySamples];
				tmpButtonLatencySeconds = [reactionTime.buttonLatencySeconds];
				
				reactionTime = [];
				reactionTime.buttonLatencySamples = tmpButtonLatencySamples;
				reactionTime.buttonLatencySeconds = tmpButtonLatencySeconds;
                    
				eventName{i}  = event(i).value;
                begsample     = event(i+2).sample - begTime*hdr.Fs;
                endsample     = event(i+2).sample + stopTime*hdr.Fs - 1;
                offset        = -begTime*hdr.Fs;
                trl(end+1, :) = round([begsample endsample offset]);
             else
                 disp('DIN1 does not follow condition trigger');
            end
        end

    end
end

[pathstr, name, ext] = fileparts(cfg.dataset);
clear name; clear ext;
outputDir  = [pathstr(1:end-3), '_dev_'];
file = fullfile(outputDir, 'reactionTimes.mat');
save(file, 'reactionTime');

results = regexp(eventName, 'c\d\d\d');
resIdx = ~cellfun('isempty', results);


eventName = eventName(resIdx);
