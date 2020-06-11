iChan = 2;
fsample = 510;
testData = data(iChan,:,1);
dft = dftmtx(length(testData));
%Do the fourier transform of the data.
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
amp(iChan,:) = abs(dftData);
cos(iChan,:) = real(dftData);
sin(iChan,:) = -imag(dftData);
freq = (fsample/2)*linspace(0,1,nfr);% freq values.  
figure;bar(freq,amp(iChan,:));

testFFT = fft(testData);
ampFFT = abs(testFFT);
ampFFTnorm = 2*ampFFT(:,1:length(ampFFT)/2+1)/length(ampFFT);
figure;bar(freq,ampFFTnorm)

% have to be done for each channel
fullData = data(:,:,1);
ftestFFT = fft(fullData);
fampFFT = abs(ftestFFT);
fampFFTnorm = 2*fampFFT(:,1:length(fampFFT)/2+1)/length(fampFFT);
figure;bar(freq,fampFFTnorm)