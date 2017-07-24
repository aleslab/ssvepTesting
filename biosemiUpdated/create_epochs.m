function [ data ] = create_epochs( cfg, datain )
%create_epochs Creates epochs. 
%   Detailed explanation goes here


% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end


cfg.trimdur    = ft_getopt(cfg, 'trimdur',    1); %Time to trim from front and back
cfg.epochdur   = ft_getopt(cfg, 'epochdur',   2); %~2s

% check if the input data is valid for this function
datain = ft_checkdata(datain, 'datatype', {'raw', 'raw+comp'}, 'hassampleinfo', 'yes');


%Lets figure out how many cycles to trim and how many to combine for an
%epoch.

%First get the size in samples of all epochs, using anonymous function and
%cell fun to get the 2 dimension of each trial.
cycleLengthSamp = cellfun(@(x) size(x,2),datain.trial);
cycleLengthSamp = unique(cycleLengthSamp);
cycleLengthSecs  = cycleLengthSamp/datain.fsample;

if length(cycleLengthSamp)>1
    warning('Cycle lengths are not equal, use resample_steadystate()');
    
    %TODO: incorporate monitor refresh info. 
    if max(diff(cycleLengthSecs))>.01667, %Why 16ms well that means cycles differ by more than a frame at 60Hz,
        error('Cycle lengths differ by more than 16 ms, cannot continue')
    end
    cycleLengthSecs = mean(cycleLengthSecs);
    cycleLengthSamp = mean(cycleLengthSamp);
end

nCyclesToTrim = round(cfg.trimdur/cycleLengthSecs); %number of cycles to remove from start and end. 
nCyclesPerEpoch = round(cfg.epochdur/cycleLengthSecs); %Number of cycles to stick together to make an "epoch" for frequency domain

trialList = unique(datain.trialinfo(:,2)); %These are are real trials

trlIdx = 1; %Lazy increment in loop FIX THIS
for iTrial = 1:length(trialList),
    thisTrial = trialList(iTrial);
    
    cyclesInThisTrial = find(datain.trialinfo(:,2)==thisTrial);
    
    %Trim the edge cycles
    trimCycleInd = cyclesInThisTrial((nCyclesToTrim+1):(end-nCyclesToTrim));
    
    nEpochs = round(length(trimCycleInd)/nCyclesPerEpoch)
    %Step through the cylces, epoch by epoch.
    for iEpoch = 1:nEpochs
        
        cycleIndices = ((iEpoch-1)*nCyclesPerEpoch + 1):iEpoch*nCyclesPerEpoch; %Count 1:n, n+1:2n, 2n+1:3n:
        cycleIndices = trimCycleInd(cycleIndices);
        data.trial{trlIdx} = cat(2,datain.trial{cycleIndices});
        data.time{trlIdx}  = (0:(length(data.trial{trlIdx})-1))./datain.fsample;
        data.trialinfo(trlIdx,1) = datain.trialinfo(cycleIndices(1),1);
        data.trialinfo(trlIdx,2) = thisTrial;
        data.trialinfo(trlIdx,3) = cycleLengthSamp;
        trlIdx = trlIdx +1;
        
    end
    
end

data.hdr = datain.hdr;
data.fsample = datain.fsample;
data.label   = datain.label;
data.cfg     = datain.cfg;

end

