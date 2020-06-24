%%% simulation of EEG response across harmonics for different duty cycles
clearvars; close all;

%%%%%%%%%%%%%%%%
%% parameters 
nT=192; 
fsample = 510;
wtSize = 1/(85/32)*1000; % in sec
% t=linspace(0,wtSize,nT); % time
sampling = wtSize/nT; % in sec
t = 0:sampling:wtSize-sampling; % time - have to remove the last point since it is the same as the first of the cycle

allDC = [0.125 0.25 0.5 0.75 0.875];
allfq = [85/8 85/16 85/32];
allTimeWin = 1./allfq*1000;
repCycle = 2;

nfr = floor(length(t)*repCycle/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
freqs = (fsample/2)*linspace(0,1,nfr);% freq values
% can also do: freqSpectrum = fsample*(0:nT)/nT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using an impulse function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impulse = (t==0|t==t(nT*allDC(dc))); % On at t0 and Off at duty-cycle

%%%%%%%%%%%%%%%%
%%%% create the response for all frequencies and duty-cycles
impWave = zeros(length(t)*repCycle,length(allfq),length(allDC));

for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        clear twin impTmp
        % create the impulse ON/OFF for a time window of the size of 1
        % frequency cycle
        twin = 0:sampling:allTimeWin(fq)-sampling; % time window for one cycle
        impTmp = (twin==0|twin==twin(length(twin)*allDC(dc)+1)); 
        % then copy it multiple times to get to the 1 sec window
        impWave(:,fq,dc) = repmat(impTmp,1,nT/length(twin)*repCycle);
    end
end


%%%% compute Fourier
impFFT = zeros(length(freqs),length(allfq),length(allDC));
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        tmpFFT = fft(impWave(:,fq,dc));
        % normalise the Fourier
        impFFT(:,fq,dc) = 2*tmpFFT(1:length(tmpFFT)/2+1)/length(tmpFFT);
    end
end


%%%% plot the waveform and Fourier spectrum
maxFreq = find(freqs==85/8*5);
for fq = 1:length(allfq)
    figure('Position', [0 0 1500 500]); hold on;
    for dc = 1 :length(allDC)
        % waveform
        subplot(2,5,dc)
%         plot(t,impulse(:,fq,dc));
        plot(impWave(:,fq,dc));
        xlabel('Time (s)')
        ylabel('Amplitude');
        title([num2str(allfq(fq)) 'Hz ' num2str(allDC(dc)*100)])
        % spectrum
        subplot(2,5,dc+5)
        bar(freqs(2:maxFreq),abs(impFFT(2:maxFreq,dc+(fq-1)*5)));
%         bar(abs(impFFT(:,dc+(fq-1)*5)));
        xlabel('Frequency (Hertz)')
        ylabel('Amplitude');
    end
%     saveas(gcf,['figures' filesep 'impulse' num2str(allfq(fq),'%.0f')],'png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using square wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% create the response for all frequencies and duty-cycles
squareWave = zeros(length(t)*repCycle,length(allfq),length(allDC));

for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        clear twin sqTmp
        twin = 0:sampling:allTimeWin(fq)-sampling; % time window for one cycle (in ms so needs to be in sec for the square function)
        sqTmp = square(2*pi*allfq(fq)*(twin/1000),allDC(dc)*100); 
        % then copy it multiple times to get to the 1 sec window
        squareWave(:,fq,dc) = repmat(sqTmp,1,nT/length(twin)*repCycle);
    end
end


%%%% compute Fourier
squareFFT = zeros(length(freqs),length(allfq),length(allDC));
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        clear tmpFFT
        tmpFFT = fft(squareWave(:,fq,dc));
        % normalise the Fourier
        squareFFT(:,fq,dc) = 2*tmpFFT(1:length(tmpFFT)/2+1)/length(tmpFFT);
    end
end


%%%% plot the waveform and Fourier spectrum
maxFreq = find(freqs==85/8*5);
for fq = 1:length(allfq)
    figure('Position', [0 0 1500 500]); hold on;
    for dc = 1 :length(allDC)
        % waveform
        subplot(2,5,dc)
%         plot(t,impulse(:,fq,dc));
        plot(squareWave(:,fq,dc));
        xlabel('Time (s)')
        ylabel('Amplitude');
        title([num2str(allfq(fq)) 'Hz ' num2str(allDC(dc)*100)])
        % spectrum
        subplot(2,5,dc+5)
        bar(freqs(2:maxFreq),abs(squareFFT(2:maxFreq,dc+(fq-1)*5)));
%         bar(abs(impFFT(:,dc+(fq-1)*5)));
        xlabel('Frequency (Hertz)')
        ylabel('Amplitude');
    end
%     saveas(gcf,['figures' filesep 'square' num2str(allfq(fq))],'png')
end




% %%%% Other way to compute FFT
% for fq = 1:length(allfq)
%     for dc = 1 :length(allDC)
%         %Do the fourier transform of the data.
%         testData = squeeze(squareWave(:,fq,dc));
%         dft = dftmtx(length(testData));
%         dftData = testData'*dft;
%         %Select just the unique frequencies.
%         nfr = floor(length(testData)/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
%         dftData = dftData(:,1:nfr);
%         freqsToDouble = 2:(nfr-1+mod(length(testData),2));
%         %Now we double the the amplitude of the freqs that are represented
%         %twice in the fourer transform
%         dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
%         %Now we normalize by the number of points in the fourier transform to
%         %get the expected amplitude value instead of the raw dot product.
%         dftData = dftData/length(testData);
%         amp(:,fq,dc) = abs(dftData);
%         %  cos = real(dftData);
%         %  sin = -imag(dftData);
%         % freq = (fsample/2)*linspace(0,1,nfr);% freq values.
%         % figure;bar(freq,amp);
%     end
% end
% 
% %%%% plot the waveform and Fourier spectrum
% maxFreq = find(freqs==85/8*5);
% for fq = 1:length(allfq)
%     figure('Position', [0 0 1500 500]); hold on;
%     for dc = 1 :length(allDC)
%         % waveform
%         subplot(2,5,dc)
% %         plot(t,impulse(:,fq,dc));
%         plot(squareWave(:,fq,dc));
%         xlabel('Time (s)')
%         ylabel('Amplitude');
%         title([num2str(allfq(fq)) 'Hz ' num2str(allDC(dc)*100)])
%         % spectrum
%         subplot(2,5,dc+5)
%         bar(freqs(1:maxFreq),abs(amp(1:maxFreq,dc+(fq-1)*5)));
% %         bar(abs(impFFT(:,dc+(fq-1)*5)));
%         xlabel('Frequency (Hertz)')
%         ylabel('Amplitude');
%     end
% %     saveas(gcf,['figures' filesep 'square' num2str(allfq(fq))],'png')
% end
