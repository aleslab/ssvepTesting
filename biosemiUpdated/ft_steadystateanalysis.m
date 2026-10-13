function [steadystate] = ft_steadystateanalysis(cfg, data)
%function [steadystate] = ft_steadystateanalysis(cfg, data)
% ft_steadystateanalysis Computes steady-state evoked potential data
%
% Use as
%   [steadystate] = ft_steadystateanalysis(cfg, data)
%
% The data should be organized so each fieldtrip "trial" is actual a
% steady-state "epoch". Each epoch should be divisable into an integer
% number of "cycles" each with identical number of samples. 
%
%
% Important note about magnitude of frequency domain coefficients: The Sine
% and Cosine components are not the coefficients of the Fourier transform.
% They are the coefficients for the Inverse fourier transform. That allows
% the coefficients to be more intuitively understood because they directly
% relate to the time domain representation of the data. 
%
%   cfg.keeptrials         = 'yes' or 'no', return individual trials or average (default = 'no')
%
% The output has the following fields:
%            
%    Data fields:
%            wave: Default cycle waveform (combined LCM cycle for two tags)
%       waveforms: Cell array of chan_time waveforms, default first
%      wavelabels: Matching GUI labels: nF1+nF2, nF1, nF2 (or nF1)
%       wavetimes: Matching time vectors in milliseconds, endpoint excluded
% wavefundamentals: Reciprocal cycle durations in Hz
% wavecyclesamples: Integer sample counts of the direct cycle averages
%        timeunit: 'ms' for time and wavetimes
% cfg.f1hz, cfg.f2hz: Tag frequencies in Hz. f2hz is optional; dual-tag epochs
%                   must contain integer cycles of both tags.
%             amp: Frequency domain amplitude.
%             sin: Frequency domain amplitude of the sin component
%             cos: Frequency domain amplitude of the cos component
%      
%             pval: The pvalue for the frequency component as calculated by the tcirc algorithm
%     stderrradius: [32x649 double]
%     confradius: The radius of the confidence interval assuming symettric error (defaults to 95%)
%             
%    Various values describing the data. 
%            i1f1: Frequency index of the first driven frequency
%           nchan: Number of channels
%         fsample: Sampling rate of time domain in Hz. 
%             nfr: Number of frequency components in 
%            ndft: Number of coefficients in fourier transform (usually (nfr-1)*2)
%              nt: Number of timepoints in time domain representation.
%            dtms: Time between samples in milliseconds.
%            dfhz: Frequency difference between spectrum values in Hz.           
%            freq: Frequency values
%            time: time values
%           label: channel labels
%      wavedimord: dimension order for wave
%          dimord: dimension order for all other data fields
%             cfg: fieldtrip config structure
%             
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end


% set the defaults
cfg.trimdur    = ft_getopt(cfg, 'trimdur',    1); %Time to trim from front and back
cfg.epochdur   = ft_getopt(cfg, 'epochdur',   2); %~2s
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials',  'no');
cfg.channel    = ft_getopt(cfg, 'channel',     'all');
cfg.f1hz       = ft_getopt(cfg, 'f1hz',  []);
cfg.f2hz       = ft_getopt(cfg, 'f2hz',  []);

% check if the input data is valid for this function
%data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp'}, 'hassampleinfo', 'yes');


% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% some proper error handling
if isfield(data, 'trial') && numel(data.trial)==0
  error('no trials were selected'); % this does not apply for MVAR data
end

if numel(data.label)==0
  error('no channels were selected');
end



%Lets figure out how many cycles to trim and how many to combine for an
%epoch.

%First get the size in samples of all epochs, using anonymous function and
%cell fun to get the 2 dimension of each trial.
%Here we're using field trip "trials" to hold each epoch. 
epochLengthSamp = cellfun(@(x) size(x,2),data.trial);
epochLengthSamp = unique(epochLengthSamp);
epochLengthSecs  = epochLengthSamp/data.fsample;

cycleLengthSamp = data.trialinfo(:,3);%TODO: Move this from trialinfo.
cycleLengthSamp = unique(cycleLengthSamp);

%TODO: Make this message more informative. 
%Need identical length trials. 
if length(epochLengthSamp)>1 || length(cycleLengthSamp)>1
    error('Epoch and/or cycle lengths are not equal, possibly use resample_steadystate() to resample trials to same length');
    
end


% The experiment design ensures complete, sample-aligned cycles.
waveCycleSamples = cycleLengthSamp;
waveLabels = {'nF1'};
if ~isempty(cfg.f2hz)
    tagCycles = round([cfg.f1hz cfg.f2hz]*epochLengthSecs);
    tagSamples = epochLengthSamp./tagCycles;
    waveCycleSamples = [lcm(tagSamples(1),tagSamples(2)) tagSamples];
    waveLabels = {'nF1+nF2','nF1','nF2'};
end
waveFundamentals = data.fsample./waveCycleSamples;

nTrials = length(data.trial);
nchan  = size(data.trial{1},1);

if (strcmp(cfg.keeptrials,'yes'))
  singleepochfourier = nan(nTrials, nchan, floor(epochLengthSamp/2)+1);
  singleepochwave = nan(nTrials, nchan, epochLengthSamp);
end


%Going to analyze each channel.
for iChan = 1:nchan,

    
    %On the first iteration let's initialize things    
    if iChan ==1
       
        steadystate.nchan = nchan;
        steadystate.fsample = data.fsample;
        steadystate.nfr = floor(epochLengthSamp/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.  
        steadystate.ndft = epochLengthSamp; %This keeps track of how long the original dft was. 
        dft = dftmtx(steadystate.ndft);        
        %dft = dft(:,1:steadystate.nfr); %Just keep the portion of the transform
        %that is unique for real signals. 
        %
        %Using a dft to analyze the data. Using dft insted of FFT for
        %historical reasons. In ancient history fft required data to be
        %power of 2 length, and could do odd things if it wasn't 
        
        %Setup the time domain representation, this is a single CYCLE, not
        %a single Epoch. 
        steadystate.nt = waveCycleSamples(1);
        dTSec = 1/data.fsample;
        steadystate.dtms  = dTSec*1000;
        steadystate.time = 0:steadystate.dtms:(steadystate.nt-1)*steadystate.dtms;
        steadystate.waveforms = cell(size(waveCycleSamples));
        steadystate.wavetimes = cell(size(waveCycleSamples));
        for iWave = 1:numel(waveCycleSamples)
            steadystate.waveforms{iWave} = NaN(nchan,waveCycleSamples(iWave));
            steadystate.wavetimes{iWave} = (0:waveCycleSamples(iWave)-1)*steadystate.dtms;
        end

        
        
        steadystate.amp  = NaN(nchan,steadystate.nfr);
        steadystate.sin  = NaN(nchan,steadystate.nfr);
        steadystate.cos  = NaN(nchan,steadystate.nfr);        
        steadystate.freq = (0:steadystate.nfr-1)*data.fsample/steadystate.ndft; % Frequency bins, including odd-length epochs.
        steadystate.dfhz = mean(diff(steadystate.freq));
        
        %TODO: Fix i1F1 spec. This just makes it up based on cycleLength or
        %finds nearest bin. 
        if isempty(cfg.f1hz)
            steadystate.i1f1 = (steadystate.ndft/cycleLengthSamp) + 1;
        else
            [~,steadystate.i1f1] = min(abs(cfg.f1hz-steadystate.freq).^2);
        end
        
        steadystate.pval   = NaN(nchan, steadystate.nfr);
        steadystate.stderrradius = NaN(nchan, steadystate.nfr);
        steadystate.confradius = NaN(nchan, steadystate.nfr);
    
    end
    
    %This is not straight forward to extract all trials and a single
    %channel from the trial cell array.  Could use a loop, but using
    %cellfun to do the work of the loop. 
    selRow = @(x) x(iChan,:); %Function to take one row of input
    thisChanData=cellfun(selRow,data.trial,'uniformoutput',false);
    thisChanData = cat(1,thisChanData{:});

    if (strcmp(cfg.keeptrials,'yes'))
        singleepochwave(:,iChan,:) = thisChanData;
    end
    
    %Do the fourier transform of the data. 
    dftData = thisChanData*dft;
    
    %Select just the unique frequencies.
    dftData = dftData(:,1:steadystate.nfr);
    
    %Now these lines are inscrutable.
    %1) The fourier transform is complex variance perserving transform
    %   That results in the coefficients not being what people expect. The
    %   frequency components that only appear once (DC, and nyquist limit)
    %   have double the value of the other coefficients
    %2) The nyquist frequency only appears when the data has
    %   an EVEN number of samples
    %       
    %The  modulus operator ise used to exclude the nyquist freq(last freq) from 
    %doubling IF an even number of samples is given. That is we go to
    %steadystate.nfr-1 if the sample is even, but to steadystate.nfr if odd. 
    freqsToDouble = 2:(steadystate.nfr-1+mod(epochLengthSamp,2)); 
    
    %Now we double the the amplitude of the freqs that are represented
    %twice in the fourer transform
    dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
    %Now we normalize by the number of points in the fourier transform to
    %get the expected amplitude value instead of the raw dot product. 
    dftData = dftData/steadystate.ndft;
   
    %For storing the single trial just keep the fourier coefficients
    %up to the user to calculate mag/sin/cos
    if (strcmp(cfg.keeptrials,'yes'))
        singleepochfourier(:,iChan,:) = dftData;
    end
    
    %Take the mean over trials of the complex valued fourier transform data
    %Note: This is IDENTICAL to taking the time down average and then
    %calculating the fourier transform.  We do the fourier first so we have
    %the single epoch coeficients we can feed into tCirc for calculating
    %statistics on the components. 
    meanDftData = mean(dftData,1);
    
    %Choice of (-) for the sin/imag component is so that in phasor diagrams
    %the convention is counter-clockwise means INCREASING latency. This is
    %identical to taking the complex conjugate of the data. But done
    %explictitly for clarity. NB: Always take care in interpreting
    %absolute phase several conventions exist in the field. 
    steadystate.amp(iChan,:) = abs(meanDftData);
    steadystate.cos(iChan,:) = real(meanDftData);
    steadystate.sin(iChan,:) = -imag(meanDftData);
        
    % Independent time-domain chain: average epochs, then average complete
    % cycles exactly as in the original single-tag implementation. No DFT
    % coefficients are used to construct any unfiltered waveform.
    meanEpoch = mean(thisChanData,1);
    for iWave = 1:numel(waveCycleSamples)
        aveWave = reshape(meanEpoch,waveCycleSamples(iWave),[])';
        aveWave = mean(aveWave,1);
        steadystate.waveforms{iWave}(iChan,:) = aveWave-mean(aveWave);
    end

    %Now lets calculate statistcs
    
    %Going to loop over the frequencies to calc the pvalue for each one 
    %The DC component is not a complex phasor so is treated differently.
    %Using a simple t-test instead of a T2 circ test.
    %
    %Removed loop by updated t2circ to accept matrix inputs
    %This results in a ~300x speedup.  The loop was causing a huge slowdown
    %When calculating many frequency values. 
    %
    %     steadystate.pval(iChan,1) = 1;
    %     steadystate.stderrradius(1) = std(dftData(:,1))/sqrt(nTrials);
    
    %Calculate DC using standard t-test. 
    
    if isreal(dftData)
        iFr = 1:steadystate.nfr; %Index into DC component
        
        [H,P,CI] = ttest(dftData(:,iFr));
        steadystate.pval(iChan,iFr) = P;
        steadystate.stderrradius(iChan,iFr) = std(dftData(:,iFr))/sqrt(nTrials);
        steadystate.confradius(iChan,iFr)   = CI(2)-mean(dftData(:,iFr));
    else
        
        iFr = 1; %Index into DC component
        
        [H,P,CI] = ttest(dftData(:,iFr));
        steadystate.pval(iChan,iFr) = P;
        steadystate.stderrradius(iChan,iFr) = std(dftData(:,iFr))/sqrt(nTrials);
        steadystate.confradius(iChan,iFr)   = CI(2)-mean(dftData(:,iFr));
        
        iFr = 2:steadystate.nfr;
        [steadystate.pval(iChan,iFr), pooledStdDev, steadystate.confradius(iChan,iFr)] = ...
            t2circ(dftData(:,iFr));
        steadystate.stderrradius(iChan,iFr) = pooledStdDev/sqrt(size(dftData,1));
    end
    
       
    
    
end

% The first direct cycle average is the backwards-compatible default.
steadystate.wavecyclesamples = waveCycleSamples;
steadystate.wavelabels = waveLabels;
steadystate.wavefundamentals = waveFundamentals;
steadystate.wave = steadystate.waveforms{1};
steadystate.time = steadystate.wavetimes{1};
steadystate.nt = size(steadystate.wave,2);
steadystate.dtms = 1000/data.fsample;
steadystate.timeunit = 'ms';
steadystate.f1hz = steadystate.freq(steadystate.i1f1);
steadystate.f2hz = cfg.f2hz;
steadystate.i1f2 = [];
if ~isempty(cfg.f2hz)
    steadystate.i1f1 = tagCycles(1)+1;
    steadystate.i1f2 = tagCycles(2)+1;
    steadystate.f1hz = steadystate.freq(steadystate.i1f1);
    steadystate.f2hz = steadystate.freq(steadystate.i1f2);
end

%If singletrial "epoch" requested
if (strcmp(cfg.keeptrials,'yes'))
    steadystate.singleepochfourier = singleepochfourier;
    steadystate.singleepochwave    = singleepochwave;
end


%allTrialMtx =cat(3,data.trial{:});
% set output variables

steadystate.label  = data.label;
steadystate.dimord = 'chan_freq'; %Applies to all fields that don't specify dimord. 
steadystate.wavedimord = 'chan_time'; %applies to waveform. 


% some fields from the input should always be copied over in the output
steadystate = copyfields(data, steadystate, {'grad', 'elec', 'opto', 'topo', 'topolabel', 'unmixing'});

if isfield(data, 'trialinfo') && strcmp(cfg.keeptrials, 'yes')
  % copy the trialinfo into the output, but not the sampleinfo
  steadystate.trialinfo = data.trialinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance steadystate
ft_postamble history    steadystate
ft_postamble savevar    steadystate
