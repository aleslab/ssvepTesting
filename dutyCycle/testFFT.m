clearvars
% fft has to be done on each channel separatly

Fs = 1000; % sampling freq
duration = 3; % seconds
timepoints = 0:1/Fs:duration-1/Fs;
freqSpectrum = Fs*(0:length(timepoints))/length(timepoints);

% create sine wave
freq1 = 5; % Hz
ampl = 2;
testData = ampl*sin(2*pi*timepoints*freq1);
figure; plot(timepoints, testData);


%%%%% Justin's method
%Do the fourier transform of the data.
dft = dftmtx(length(testData));
dftData = testData*dft;
%Select just the unique frequencies.
nfr = floor(length(testData)/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.  
dftData = dftData(:,1:nfr);
freqsToDouble = 2:(nfr-1+mod(length(testData),2));
%Now we double the the amplitude of the freqs that are represented
%twice in the fourer transform
dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
%Now we normalize by the number of points in the fourier transform to
%get the expected amplitude value instead of the raw dot product.
dftData = dftData/length(testData);
amp = abs(dftData);
cos = real(dftData);
sin = -imag(dftData);
freq = (Fs/2)*linspace(0,1,nfr);% freq values.  
figure;bar(freq,amp);


%%%%% Rama's method
otherFFT = fft(testData);
ampFFT = abs(otherFFT);
ampFFTnorm = 2*ampFFT(1:length(ampFFT)/2+1)/length(ampFFT);
figure;bar(freq,ampFFTnorm)

% go back to testData from fft
backOther = ifft(otherFFT);
test =dftData/dft(:,1:nfr);

back = ifft(complex(cos,sin));
ampl * sin(timepoints + angle(testData))