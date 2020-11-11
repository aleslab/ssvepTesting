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

allDC = [0.125 0.25 0.5 0.75 0.875]; % high res: allDC = [0.01:0.01:0.99];
allfq = [85/8 85/16 85/32];
allTimeWin = 1./allfq*1000;
repCycle = 6; % this corresponds to the data

nfr = floor(length(t)*repCycle/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
freqs = (fsample/2)*linspace(0,1,nfr);% freq values
% can also do: freqSpectrum = fsample*(0:nT)/nT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using an impulse function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impulse = (t==0|t==t(nT*allDC(dc))); % On at t0 and Off at duty-cycle

%%%%%%%%%%%%%%%%
%%%% could do with +1 Onset and -1 offset to get back the other harmonics
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
        plot(impWave(:,fq,dc));
        xlabel('Time (ms)')
        ylabel('Amplitude');
        title([num2str(allfq(fq)) 'Hz ' num2str(allDC(dc)*100)])
        % spectrum
        subplot(2,5,dc+5)
        bar(freqs(2:maxFreq),squeeze(abs(impFFT(2:maxFreq,fq,dc))));
        xlabel('Frequency (Hertz)')
        ylabel('Amplitude');
    end
    saveas(gcf,['figures' filesep 'impulse' num2str(allfq(fq),'%.0f')],'png')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using square wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% create the response for all frequencies and duty-cycles
squareWave = zeros(length(t)*repCycle,length(allfq),length(allDC));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ATTENTION!!!!!!!!!!!!!!!!!!!!!
% % square function is not doing a proper work for 5Hz! 
% % use repelem instead
% for fq = 1:length(allfq)
%     for dc = 1 :length(allDC)
%         clear twin sqTmp
%         twin = 0:sampling:allTimeWin(fq)-sampling; % time window for one cycle (in ms so needs to be in sec for the square function)
%         sqTmp = square(2*pi*allfq(fq)*(twin/1000),allDC(dc)*100); 
%         figure;plot(sqTmp)
%         % then copy it multiple times to get to the 1 sec window
%         squareWave(:,fq,dc) = repmat(sqTmp,1,nT/length(twin)*repCycle);
%     end
% end

% number of time points per frequency cycle
for fq = 1:length(allfq)
    sizeCycle(fq) = length(0:sampling:allTimeWin(fq)-sampling);
end
a = [1 -1];
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        sqTmp = repelem(a,[sizeCycle(fq)*allDC(dc) sizeCycle(fq)*(1-allDC(dc))]);
%         figure;plot(sqTmp)
        % then copy it multiple times to get to the 1 sec window
        squareWave(:,fq,dc) = repmat(sqTmp,1,nT/sizeCycle(fq)*repCycle);
    end
end
% % round repelem for high res
% for fq = 1:length(allfq)
%     for dc = 1 :length(allDC)
%         sqTmp = repelem(a,[round(sizeCycle(fq)*allDC(dc)) round(sizeCycle(fq)*(1-allDC(dc)))]);
% %         figure;plot(sqTmp)
%         % then copy it multiple times to get to the 1 sec window
%         squareWave(:,fq,dc) = repmat(sqTmp,1,nT/sizeCycle(fq)*repCycle);
%     end
% end

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
        plot(squareWave(:,fq,dc));
        xlabel('Time (ms)')
        ylabel('Amplitude');
        title([num2str(allfq(fq)) 'Hz ' num2str(allDC(dc)*100)])
        % spectrum
        subplot(2,5,dc+5)
        bar(freqs(2:maxFreq),squeeze(abs(squareFFT(2:maxFreq,fq,dc))));
        xlabel('Frequency (Hertz)')
        ylabel('Amplitude');
    end
    saveas(gcf,['figures' filesep 'square' num2str(allfq(fq))],'png')
end



%%%%% for the expt I have to add the motion conditions and reshape to get
%%%%% 22 conditions
% last conditions are (fq DC): 2.5 12.5; 5 25; 10 50; 5 75; 2.5 87.5; 2.5
% 50; 5 50;
% ATTENTION! fq1 = 10Hz fq3 = 2.5Hz
sqTmp = permute(squareWave,[1 3 2]);
sqAll = reshape(sqTmp,[1152,15]); % high res: sqAll = reshape(sqTmp,[1152,99*3]);
sqAll(:,16:22) = [squareWave(:,3,1) squareWave(:,2,2) squareWave(:,1,3) ...
    squareWave(:,2,4) squareWave(:,3,5) squareWave(:,3,3) squareWave(:,2,3)];

%%%% compute Fourier
sqAllFFT = zeros(length(freqs),size(sqAll,2));
for cond = 1:size(sqAllFFT,2)
    clear tmpFFT
    tmpFFT = fft(sqAll(:,cond));
    % normalise the Fourier
    sqAllFFT(:,cond) = abs(2*tmpFFT(1:length(tmpFFT)/2+1)/length(tmpFFT));
end
save('squareFFT.mat','sqAllFFT')


%%% plot as in the experiment: sum harmonics for the 5 DC
figure;
subplot(2,1,1);hold on
harm1 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(1)));
plot(allDC,sum(sqAllFFT(harm1,1:5),1),'b','LineWidth',2)
harm2 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(2)));
plot(allDC,sum(sqAllFFT(harm2,6:10),1),'r','LineWidth',2)
harm3 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(3)));
plot(allDC,sum(sqAllFFT(harm3,11:15),1),'g','LineWidth',2)
xlim([0 1])
legend('10Hz','5Hz','2Hz')
xlabel('duty cycle')
ylabel('sum fq amp < 50Hz')
title('sum of amplitudes')
subplot(2,1,2);hold on
plot(allDC,sum(sqAllFFT(harm1,1:5).^2),'b','LineWidth',2)
plot(allDC,sum(sqAllFFT(harm2,6:10).^2),'r','LineWidth',2)
plot(allDC,sum(sqAllFFT(harm3,11:15).^2,1),'g','LineWidth',2)
title('sum of power (squared amplitudes)')

saveas(gcf,['figures' filesep 'normalisation'],'png')



% %%%% plot high resolution
% squareFFT(~abs(squareFFT)) = NaN;
% maxFreq = find(freqs==85/8*5);
% toplot = 42:2:56;
% for fq = 1:length(allfq)
%     figure('Position', [0 0 1500 1500]); hold on;
%     for dc = 1:length(toplot)
%         % spectrum
%         subplot(4,2,dc)
% %         bar(freqs(2:maxFreq),squeeze(abs(squareFFT(2:maxFreq,fq,toplot(dc)))));
%         stem(freqs(2:end),squeeze(abs(squareFFT(2:end,fq,toplot(dc)))));
%         xlim([0 freqs(end)])
%         xlabel('Frequency (Hertz)')
%         ylabel('Amplitude');
%         title([num2str(allfq(fq)) 'Hz ' num2str(toplot(dc)) '%'])
%     end
%     saveas(gcf,['figures' filesep 'FourierAll' num2str(round(allfq(fq)))],'png')
% end
% 
% 
% figure; hold on
% harm1 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(1)));
% plot(allDC,sum(sqAllFFT(harm1,1:99),1),'b','LineWidth',2)
% harm2 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(2)));
% plot(allDC,sum(sqAllFFT(harm2,100:198),1),'r','LineWidth',2)
% harm3 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(3)));
% plot(allDC,sum(sqAllFFT(harm3,199:297),1),'g','LineWidth',2)
% xlim([0 1])
% legend('10Hz','5Hz','2Hz')
% xlabel('duty cycle')
% ylabel('sum fq amp < 50Hz')
% title('normalisation')
% saveas(gcf,['figures' filesep 'normalisationHighRes'],'png')






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
