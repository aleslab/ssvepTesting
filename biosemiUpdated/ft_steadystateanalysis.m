function [steadystate] = ft_steadystateanalysis(cfg, data)


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


% set the defaults
cfg.trimdur    = ft_getopt(cfg, 'trimdur',    1); %Time to trim from front and back
cfg.epochdur   = ft_getopt(cfg, 'epochdur',   2); %~2s
cfg.keeptrials   = ft_getopt(cfg, 'keeptrials',  'no');
cfg.channel      = ft_getopt(cfg, 'channel',     'all');

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


nTrials = length(data.trial);
nChan  = size(data.trial{1},1);


%Going to analyze each channel.
for iChan = 1:nChan,

    
    %On the first iteration let's initialize things    
    if iChan ==1
       
        Axx.fsample = data.fsample;
        Axx.nFr = floor(epochLengthSamp/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.  
        dft = dftmtx(epochLengthSamp); 
        %Using a dft to analyze the data. Using dft insted of FFT for
        %historical reasons. In ancient history fft required data to be
        %power of 2 length, and could do odd things if it wasn't 
        
        %Setup the time domain representation, this is a single CYCLE, not
        %a singly Epoch. 
        Axx.nT = cycleLengthSamp;
        Axx.dTSec = 1/data.fsample;
        Axx.dTms  = Axx.dTSec*1000;
        Axx.time = 0:Axx.dTSec:(Axx.nT-1)*Axx.dTSec;
        Axx.Wave = NaN(nChan,cycleLengthSamp);
        
        Axx.Amp  = NaN(nChan,Axx.nFr);
        Axx.Sin  = NaN(nChan,Axx.nFr);
        Axx.Cos  = NaN(nChan,Axx.nFr);        
        Axx.freq = (data.fsample/2)*linspace(0,1,Axx.nFr);% freq values.  
        Axx.dFhz = mean(diff(Axx.freq));
        
        Axx.tCircPval   = NaN(nChan, Axx.nFr);
        Axx.tCircStdErr = NaN(nChan, Axx.nFr);
    
    end
    
    %This is not straight forward to extract all trials and a single
    %channel from the trial cell array.  Could use a loop, but using
    %cellfun to do the work of the loop. 
    selRow = @(x) x(iChan,:); %Function to take one row of input
    thisChanData=cellfun(selRow,data.trial,'uniformoutput',false);
    thisChanData = cat(1,thisChanData{:});
    
    
    %Do the fourier transform of the data. 
    dftData = thisChanData*dft;
    
    %Select just the unique frequencies.
    dftData = dftData(:,1:Axx.nFr);
    
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
    %Axx.nFr-1 if the sample is even, but to Axx.nFr if odd. 
    freqsToDouble = 2:(Axx.nFr-1+mod(epochLengthSamp,2)); 
    
    %Now we double the the amplitude of the freqs that are represented
    %twice in the fourer transform
    dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
    %Now we normalize by the number of points in the fourier transform to
    %get the expected amplitude value instead of the raw dot product. 
    dftData = dftData/Axx.nFr;
   
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
    Axx.Amp(iChan,:) = abs(meanDftData);
    Axx.Cos(iChan,:) = real(meanDftData);
    Axx.Sin(iChan,:) = -imag(meanDftData);
        
    %Let's now create the time-domain representation
    %For this we are going to make a waveform that only 1 cycle long
    %instead of 1 whole epoch long. A single cycle represents the waveform
    %as the most intuitive representation of the time-domain information
    %and variabillity. 
    
    %First average all epochs together, then average all cycles together.
    %This takes being very careful in how matrices are reshaped and
    %averaged over. The exact dimension order and transpose is important.
    aveWave=mean(thisChanData,1);
    aveWave=reshape(aveWave,cycleLengthSamp,[])';
    aveWave=mean(aveWave,1);

    %Now assign the mean cycle waveform to the time domain average:
    %We are going to remove the mean.  We almost never want to plot it. If
    %we need it it's available in the DC component of Cos. 
    Axx.Wave(iChan,:)= aveWave-mean(aveWave);
    
    %Now lets calculate statistcs
    
    %Going to loop over the frequencies to calc the pvalue for each one 
    %The DC component is not a complex phasor so is treated differently.
    
%     Axx.tCircPval(iChan,1) = 1;
%     Axx.tCircStdErr(1) = std(dftData(:,1))/sqrt(nTrials);
    for iFr = 1:Axx.nFr,
        
        if isreal(dftData(:,iFr)),
            Axx.tCircPval(iChan,iFr) = 1;
            Axx.tCircStdErr(iChan,iFr) = std(dftData(:,iFr))/sqrt(nTrials);
            continue
        end
        
        [Axx.tCircPval(iChan,iFr)  Axx.tCircStdErr(iChan,iFr)] = tcirc(dftData(:,iFr));
    end
    
    
    
    
end

%allTrialMtx =cat(3,data.trial{:});
% set output variables
steadystate        = Axx;
steadystate.label  = data.label;
steadystate.dimord = 'chan_freq';


% some fields from the input should always be copied over in the output
steadystate = copyfields(data, steadystate, {'grad', 'elec', 'opto', 'topo', 'topolabel', 'unmixing'});

if isfield(data, 'trialinfo') && strcmp(cfg.keeptrials, 'yes')
  % copy the trialinfo into the output, but not the sampleinfo
  steadystate.trialinfo = data.trialinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance steadystate
ft_postamble history    steadystate
ft_postamble savevar    steadystate
