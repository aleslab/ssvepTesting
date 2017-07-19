% clean data?
% reference 32=Cz? 16=Oz?
% 1 frame difference with the photodiode - any readjustment?
% baseline window? No need? (detrending?)
% 


clear all

% addpath D:\GitHubRepo\fieldtrip
% ft_defaults
% cfg.dataset   = 'D:\StAndrews\testSSVEP\SSVEPtest3.bdf';

addpath /Users/marleneponcet/Documents/Git/fieldtrip
ft_defaults
cfg.dataset = '/Users/marleneponcet/Documents/ssvep/test/testSSVEP/SSVEPtestLR.bdf';

%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.padding      = length (in seconds) to which the trials are padded for filtering (default = 0)
%   cfg.padtype      = string, type of padding (default: 'data' padding or
%                      'mirror', depending on feasibility)
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%

cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.trialfun = 'lock2SsvepTag_LR'; 

[cfg] = ft_definetrial(cfg);

cfg.demean        ='yes';
cfg.reref         = 'yes';
cfg.refchannel    = {'A32'}; % 32=Cz?
[data] = ft_preprocessing(cfg);

% %%% A bit of visualization
% % plot the same channel over all epochs ( = cycles + trials)
% figure; hold on;
% for nbtrials = 1:5%length(data.trial)
%     test = data.trial{1,nbtrials};
%     plot(test(1,:),'Color',[rand(1) rand(1) rand(1)])
% end

cfg = [];
[timelock] = ft_timelockanalysis(cfg, data); % computes the timelocked average ERP/ERF

% plot the average for a specified channel
% plot(timelock.avg(1,:))

%%% Power Spectrum
fftDat = fft(timelock.avg(1:32,:)'); % for the 32 channels
% spectrumDat = abs(2*fftDat(2:floor(end/2))/length(fftDat));
spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric). 
% Also transforms the matrix in a single vector. To keep the channels
% separate (2D matrix): spectrumDat = abs(fftDat(2:floor(end/2),:));
spectrumDat = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
freqs = [1:length(spectrumDat)] * 1/(timelock.time(end)-timelock.time(1));
figure; hold on; plot(freqs(freqs<80),spectrumDat(freqs<80),'r');
pdSpecPlot_copy(freqs(1:80),spectrumDat(1:80),[]);

spectrumChan = abs(fftDat(2:floor(end/2),16)); % for the specified channel, here 16=Oz?
spectrumChan = spectrumChan*2/length(spectrumChan);
figure; hold on; plot(freqs(freqs<80),spectrumChan(freqs<80),'r');
pdSpecPlot_copy(freqs(1:80),spectrumChan(1:80),[]);


%%% Below uses dfmtx BUT
%%%   If X is a column vector of length N, then
%%%   DFTMTX(N)*X yields the same result as FFT(X); however, 
%%%   FFT(X) is more efficient.
% Axx.Wave = timelock.avg(1:32,:)';
% Axx.nT = size(Axx.Wave,1);
% Axx.nFr = round(size(Axx.Wave,1)/2);
% dft = dftmtx(Axx.nT);
% 
% dftDat = dft*Axx.Wave;
% dftDat = dftDat(1:Axx.nFr,:);
% 
% Axx.dFHz = (data.hdr.Fs)/Axx.nT;
% 
% % Stuff tends to be multiples of 1 Hz. But we don't have that info
% % here so we are just going to set the index to 1 to make nF1 be all
% % freqs
% Axx.i1F1       = 1;
% Axx.i1F2       = 0;
% 
% Axx.Amp = abs(2*(dftDat/Axx.nFr));
% Axx.Cos = 2*real(dftDat)/Axx.nFr;
% Axx.Sin = -2*imag(dftDat)/Axx.nFr;
% 
% 
% 
% freqs = 0:Axx.dFHz:(Axx.dFHz*(Axx.nFr-1));
% 
% pdSpecPlot(freqs(1:80),Axx.Amp(1:80,16)',[])



%%%%%%%%%%%%%%%%%%
%%%% depending on condition (saved in data.trialinfo)
%%%%%%%%%%%%%%%%%%

allcond = unique(data.trialinfo);
for cond = 1 : length(allcond)
    cfg = [];
    cfg.trials = find(data.trialinfo == allcond(cond));
    condAverage(cond).data = ft_timelockanalysis(cfg, data);
    fftDat = fft(condAverage(cond).data.avg(16,:)'); % for the 32 channels
    spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric).
    condSpect(cond,:) = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
end

freqs = [1:length(spectrumDat)] * 1/(condAverage(cond).data.time(end)-condAverage(cond).data.time(1));

figure; hold on; 
for cond = 1 : length(allcond)
    subplot(3,3,cond)
    plot(freqs(freqs<40),condSpect(cond,freqs<40),'r');
end