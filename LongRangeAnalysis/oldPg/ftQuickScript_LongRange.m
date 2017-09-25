% still have to figure out how to remove noise and bad channels
% A23 = Oz

clear cfg data;

addpath /Users/marleneponcet/Documents/Git/fieldtrip
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults
cfg.dataset   =  '/Users/marleneponcet/Documents/LongRangeSSVEP/LongRangeS01.bdf';

% cfg.dataset   =  'D:\StAndrews\LongRangeSSVEP\LongRangeS01.bdf'
% addpath D:\GitHubRepo\fieldtrip
% addpath D:\GitHubRepo\ssvepTesting\svndlCopy
% addpath D:\GitHubRepo\ssvepTesting\biosemiUpdated

% cfg.dataset   = '/Volumes/Amrutam/Marlene/LongRangeSSVEP/LongRangeS01.bdf'
% addpath /Volumes/Amrutam/Marlene/Git/fieldtrip
% addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
% addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated

% define trials then preprocessing or the other way around?? Faster this
% way.. 

cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.channel = 1:128;
cfg.trialdef.preStimDuration = 1.2;
cfg.trialdef.nbTotalCycles = 31;
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'lock2ssvep_LongRange'; 
[cfg] = ft_definetrial(cfg);


cfg.lpfilter  = 'no';
cfg.lpfreq        = 100;
cfg.demean        ='yes';
cfg.reref         = 'no'; 
cfg.refchannel    = {'A3'}; % A3 = CPz
[data] = ft_preprocessing(cfg); 

cfg = [];
cfg.vartrllength = 2;
[timelockEpoch] = ft_timelockanalysis(cfg, data)

%remove the extra channels. 
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
dataEeg = ft_selectdata(cfg,data);
%Do the steadystate analysis
[Axx] = ft_steadystateanalysis(cfg, data)


%Now do a quick spectrum of the 12th channel:
figure;
pdSpecPlot(Axx.freq(2:80),Axx.Amp(12,2:80)',Axx.tCircPval(12,2:80)<.05)
title('Channel 12')

%Now setup an interactive plot
cfg.layout = 'biosemi128.lay';
cfg.channel = {'all'};
figure;clf;
interactiveTopoSpecPlot(cfg,Axx)
interactiveSteadyStatePlot(cfg,Axx)


% %%% Power Spectrum
% fftDat = fft(timelockEpoch.avg(1:128,:)'); % for the 32 channels
% % spectrumDat = abs(2*fftDat(2:floor(end/2))/length(fftDat));
% spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric). 
% % Also transforms the matrix in a single vector. To keep the channels
% % separate (2D matrix): spectrumDat = abs(fftDat(2:floor(end/2),:));
% spectrumDat = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
% freqs = [1:length(spectrumDat)] * 1/(timelockEpoch.time(end)-timelockEpoch.time(1));
% figure; hold on; plot(freqs(freqs<80),spectrumDat(freqs<80),'r');
% pdSpecPlot(freqs(1:80),spectrumDat(1:80),[]);
% 
% spectrumChan = abs(fftDat(2:floor(end/2),16)); % for the specified channel, here 16=Oz?
% spectrumChan = spectrumChan*2/length(spectrumChan);
% figure; hold on; plot(freqs(freqs<80),spectrumChan(freqs<80),'r');
% pdSpecPlot(freqs(1:80),spectrumChan(1:80),[]);


%%%%%%%%%%%%%%%%%%
%%%% depending on condition (saved in data.trialinfo)
%%%%%%%%%%%%%%%%%%

allcond = unique(data.trialinfo(:,1));
for cond = 1 : length(allcond)
    cfg = [];
    cfg.trials = find(data.trialinfo(:,1) == allcond(cond));
    condAverage(cond).data = ft_timelockanalysis(cfg, data);
    fftDat = fft(condAverage(cond).data.avg(23,:)'); 
    spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric).
    condSpect(cond,:) = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
end

freqs = [1:length(spectrumDat)] * 1/(condAverage(cond).data.time(end)-condAverage(cond).data.time(1));

figure; hold on; 
for cond = 1 : length(allcond)
    subplot(3,3,cond)
    plot(freqs(freqs<40),condSpect(cond,freqs<40),'r');
end



%%%%% follow up
allcond = unique(data.trialinfo(:,1));
for cond=1:length(allcond)
    cfg.trials = find(data.trialinfo(:,1) == allcond(cond));
%     cfg.channel =  {'all','-GSR1', '-GSR2', '-Erg1','-Erg2','-Resp','-Plet','-Temp', '-Status'};
    % dataEeg = ft_selectdata(cfg,data);
    [Axx(cond)] = ft_steadystateanalysis(cfg, data)
%     cfg.channel = {'all'};
    figure;clf;
    % interactiveTopoSpecPlot(cfg,Axx(cond))
    interactiveSteadyStatePlot(cfg,Axx(cond))
end