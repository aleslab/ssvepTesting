%%% simulation of a sinusoid response for different stimulus stimulations 
%%% stimulations are modelled as squarewave or impulse function
%%% This is to test the effect of frequency stimulation and duty-cycle on
%%% fundamental and harmonics SSVEP amplitudes

clearvars; close all;
% addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated/')

%%%%%%%%%%%%%%%%
%% parameters 
fsample = 2048*2; % sampling rate (number of timepoints). Needs to be high to be able to get high resolution of DC
wtSize = 2; % size of the time window in sec. Should be long enough to fit the lowest fq (here would be 1Hz)
t=linspace(0,wtSize-1/fsample,fsample); % time
sampling = wtSize/fsample; % in sec

allDC = [0.125 0.25 0.5 0.75 0.875]; % 0.1:0.1:0.9; % duty-cycles
% lowest fq limited by wtSize, highest limited by fsample & how duty-cycle
% cuts the waveform!
allfq = [1 4 8 12 16]/wtSize; % 1:20/wtSize; % frequency rate

nfr = floor(length(t)/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
freqs = (fsample/2)*linspace(0,1,nfr);% freq values 


%%%%%%%%%%%%%%%%
%% square wave
%%%% initialise squarewave for all frequencies and duty-cycles
squareWave = zeros(length(t),length(allfq),length(allDC));

% number of time points corresponding of 1 cycle (different for different
% fq rate)
for fq = 1:length(allfq)
    sizeCycle(fq) = length(0:sampling:1/allfq(fq)-sampling);
end

% create squarewave + impulse stimulation for all conditions
% also create a "trigger" stim ON used for the ERP
a = [1 0]; % squarewave limits = 1 0
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        % create one cycle
        sqTmp = repelem(a,single([sizeCycle(fq)*allDC(dc) sizeCycle(fq)*(1-allDC(dc))]));
        % repeat the cycle to fit the 1 sec window
        squareWave(:,fq,dc) = repmat(sqTmp,1,fsample/sizeCycle(fq));
        % use the squarewave to create impulse & trigger
        squareWave2 = [0 squareWave(:,fq,dc)'];
        impulseStim(:,fq,dc) = abs(squareWave2(1:end-1) - squareWave(:,fq,dc)');
        triggerON(:,fq,dc) = squareWave2(1:end-1) - squareWave(:,fq,dc)'<0;
    end
end



%%%% compute Fourier
squareFFT = zeros(length(freqs),length(allfq),length(allDC));
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        clear tmpFFT
        tmpFFT = fft(squareWave(:,fq,dc));
        % normalise the Fourier
        squareFFT(:,fq,dc) = abs(2*tmpFFT(1:length(tmpFFT)/2+1)/length(tmpFFT));
    end
end



% %%%%%%%%%%%%%%%%
%% ERP
% single ERP 
nbSin = 4; % Sinwave will be of 4 times smaller than the entire window (= 250 ms ERP)
tsmall = t(1:length(t)/nbSin); 
erp = sin(nbSin*2*pi*tsmall);
% set to full window
erpFull = [erp zeros(1,length(erp)*(nbSin-1))];

% create SSVEP for different fq rate and DC
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        timeON = find(triggerON(:,fq,dc) == 1)-1; % remove 2 for the circshift 
        % do arrayfun?
        test=arrayfun(@(x) circshift(erpFull,timeON(x)), 1:length(timeON), 'UniformOutput', false );
        erpWave(:,fq,dc) = sum(reshape(cell2mat(test),length(t),length(timeON)),2)';
%         figure;plot(t,erpWave(:,fq,dc))
    end
end

figure;plot(erp)
ftERP = fft(erp);

