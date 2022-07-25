% Compare the results of recreating freq amp from the cycle vs the entire
% epoch

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults


load('/Users/marleneponcet/Documents/data/dutyCycle/cleanData/DutyCycle_newStim_S02_20180503_clean.mat')

cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
cfg.layout = 'biosemi128.lay';

% select just 1 condition to work with
cfg.trials = find(cleanData.trialinfo(:,1) == 108); 
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
data = ft_selectdata(tmpcfg, cleanData);



%% FROM FT_STEADYSTATEANALYSIS

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

nTrials = length(data.trial);
nchan  = size(data.trial{1},1);

%Going to analyze each channel.
for iChan = 1:nchan

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
        steadystate.nt = cycleLengthSamp;
        dTSec = 1/data.fsample;
        steadystate.dtms  = dTSec*1000;
        steadystate.time = 0:steadystate.dtms:(steadystate.nt-1)*steadystate.dtms;
        steadystate.wave = NaN(nchan,cycleLengthSamp);

        steadystate.amp  = NaN(nchan,steadystate.nfr);
        steadystate.sin  = NaN(nchan,steadystate.nfr);
        steadystate.cos  = NaN(nchan,steadystate.nfr);        
        steadystate.freq = (data.fsample/2)*linspace(0,1,steadystate.nfr);% freq values.  
        steadystate.dfhz = mean(diff(steadystate.freq));
        
        steadystate.i1f1 = (steadystate.ndft/cycleLengthSamp) + 1;
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
        
    %Let's now create the time-domain representation
    %For this we are going to make a waveform that only 1 cycle long
    %instead of 1 whole epoch long. A single cycle represents the waveform
    %as the most intuitive representation of the time-domain information
    %and variabillity. 
    
    % keep the full epoch and create a mean cycle
    aveWave=mean(thisChanData,1);
    steadystate.avewave(iChan,:)= aveWave-mean(aveWave);
    
    cycleWave = reshape(aveWave,cycleLengthSamp,[])';
    cycleWave = mean(cycleWave,1);
    steadystate.cycleWave(iChan,:)= cycleWave-mean(cycleWave);
    
end







%%%%
figure;bar(steadystate.freq(2:100),steadystate.amp(19,2:100));title('original')

% first try for just a summated trial
pickChan = 19; % pick one channel so that it's faster
summatedCycle.wave = steadystate.cycleWave(pickChan,:) + steadystate.cycleWave(pickChan,:);
summatedCycle.sin = steadystate.sin(pickChan,:) + steadystate.sin(pickChan,:);
summatedCycle.cos = steadystate.cos(pickChan,:) + steadystate.cos(pickChan,:);
summatedCycle.amp = sqrt(summatedCycle.sin.^2 + summatedCycle.cos.^2);
figure;bar(steadystate.freq(2:100),summatedCycle.amp(1,2:100));title('summated')


summatedEpoch = steadystate.avewave(pickChan,:) + steadystate.avewave(pickChan,:);
dftData = summatedEpoch*dft;
dftData = dftData(:,1:steadystate.nfr);
freqsToDouble = 2:(steadystate.nfr-1+mod(epochLengthSamp,2));
dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
dftData = dftData/steadystate.ndft;
meanDftData = mean(dftData,1);
summatedEpochAmp = abs(meanDftData);
figure;bar(steadystate.freq(2:100),summatedEpochAmp(1,2:100))


% shift wave
shift = circshift(steadystate.cycleWave(pickChan,:),[0 length(steadystate.cycleWave(pickChan,:))/2]); 
cycleShift.wave = steadystate.cycleWave(pickChan,:) + shift;
figure;plot(steadystate.time,steadystate.cycleWave(pickChan,:));hold on;plot(steadystate.time,cycleShift.wave);
plot(steadystate.time,shift);

periods = 1 ./steadystate.freq;
phaseShift = steadystate.time(end/2+1)/1000 ./ periods .* 360; 
sinH = steadystate.sin(pickChan,:);
cosH = steadystate.cos(pickChan,:);
shiftSin = sind(phaseShift) .*  cosH + cosd(phaseShift) .* sinH;
shiftCos = cosd(phaseShift) .* cosH - sind(phaseShift) .* sinH ;

cycleShift.sin = sinH + shiftSin;
cycleShift.cos = cosH + shiftCos;
cycleShift.amp = sqrt(cycleShift.sin.^2 + cycleShift.cos.^2);
figure;bar(steadystate.freq(2:100),cycleShift.amp(2:100));title('cycle')

% reconstruct long signal with shifts every half cycle
for cc=1: length(steadystate.avewave) / size(steadystate.cycleWave,2) % nb of cycles in the epoch
    cycleTmp = steadystate.avewave(pickChan,size(steadystate.cycleWave,2)*(cc-1)+1 : size(steadystate.cycleWave,2)*cc);
    shiftTmp = circshift(cycleTmp,[0 length(cycleTmp)/2]); 
    recSig(size(steadystate.cycleWave,2)*(cc-1)+1:size(steadystate.cycleWave,2)*cc) = shiftTmp+cycleTmp;
%     figure;plot(steadystate.time,cycleTmp);hold on;plot(steadystate.time,recSig(size(steadystate.cycleWave,2)*(cc-1)+1:size(steadystate.cycleWave,2)*cc)); plot(steadystate.time,shiftTmp);
end
figure;plot(recSig); hold on; plot(steadystate.avewave(pickChan,:))
cycleFromEpoch = reshape(recSig,cycleLengthSamp,[])';
figure;plot(mean(cycleFromEpoch));hold on;plot(cycleShift.wave);  

clear dftData
dftData = recSig*dft;
dftData = dftData(:,1:steadystate.nfr);
freqsToDouble = 2:(steadystate.nfr-1+mod(epochLengthSamp,2));
dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
dftData = dftData/steadystate.ndft;
meanDftData = mean(dftData,1);
shiftEpochAmp = abs(meanDftData);
figure;bar(steadystate.freq(2:100),shiftEpochAmp(1,2:100));title('epoch')

figure;bar(steadystate.freq(2:100),shiftEpochAmp(1,2:100));hold on; bar(steadystate.freq(2:100),cycleShift.amp(2:100),'LineStyle','-','LineWidth',1.5,'EdgeColor','r','FaceColor','none');

