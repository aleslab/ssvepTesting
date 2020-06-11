Fs = 1000; % sampling freq
duration = 3; % seconds
timepoints = 0:1/Fs:duration-1/Fs;
freqSpectrum = Fs*(0:length(timepoints))/length(timepoints);

% create sine wave
freq1 = 5; % Hz
ampl = 2;
wav = ampl*sin(2*pi*timepoints*freq1);
figure; plot(timepoints, wav);

% compute fft
ourfft = fft(wav);
allampl = abs(ourfft);
figure; plot(freqSpectrum(1:end-1),allampl);

% plot normalised amplitudes
normampl = 2*allampl(1:length(timepoints)/2+1)/length(timepoints);
figure; stem(freqSpectrum(1:50), normampl(1:50));

% plot phase
fftForPhase = ourfft;
fftForPhase(abs(ourfft)<1e-6) = 0;
allphase = angle(fftForPhase);
figure;stem(freqSpectrum(1:50),rad2deg(allphase(1:50)));

%% modifying the sine wave
% create sine wave
freq1 = 5; % Hz
ampl1 = 1.2;
freq2 = 7; % Hz
ampl2 = 0.5;
freq3 = 11; % Hz
ampl3 = 2.5;
wav = ampl1*sin(2*pi*timepoints*freq1)+ampl2*sin(2*pi*timepoints*freq2)+ampl3*sin(2*pi*timepoints*freq3);
% wav = wav+0.5*randn(size(wav)); % add noise
figure; plot(timepoints, wav);

% compute fft
ourfft = fft(wav);
allampl = abs(ourfft);

% plot normalised amplitudes
normampl = 2*allampl(1:length(timepoints)/2+1)/length(timepoints);
figure; stem(freqSpectrum(1:50), normampl(1:50));

% plot phase
fftForPhase = ourfft;
fftForPhase(abs(ourfft)<1e-6) = 0;
allphase = angle(fftForPhase);
figure;stem(freqSpectrum(1:50),allphase(1:50)/pi)