%%% simulation of EEG response across harmonics for different duty cycles
clearvars; close all;

%%%%%%%%%%%%%%%%
%% parameters 
nT=192*1; 
fsample = 510*1;
wtSize = 1/(85/32)*1000; % in sec
% t=linspace(0,wtSize,nT); % time
sampling = wtSize/nT; % in sec
t = 0:sampling:wtSize-sampling; % time - have to remove the last point since it is the same as the first of the cycle

% low res:
% allDC = [0.125 0.25 0.5 0.75 0.875]; 
% high res: 
allDC = [0.01:0.01:0.99];
allfq = [85/8 85/16 85/32];
allTimeWin = 1./allfq*1000;
repCycle = 6; % this corresponds to the data

nfr = floor(length(t)*repCycle/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
freqs = (fsample/2)*linspace(0,1,nfr);% freq values
% can also do: freqSpectrum = fsample*(0:nT)/nT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using square wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% create the response for all frequencies and duty-cycles
squareWave = zeros(length(t)*repCycle,length(allfq),length(allDC));

% number of time points per frequency cycle
for fq = 1:length(allfq)
    sizeCycle(fq) = length(0:sampling:allTimeWin(fq)-sampling);
end
a = [1 -1];
% for fq = 1:length(allfq)
%     for dc = 1 :length(allDC)
%         sqTmp = repelem(a,[sizeCycle(fq)*allDC(dc) sizeCycle(fq)*(1-allDC(dc))]);
% %         figure;plot(sqTmp)
%         % then copy it multiple times to get to the 1 sec window
%         squareWave(:,fq,dc) = repmat(sqTmp,1,nT/sizeCycle(fq)*repCycle);
%     end
% end
% round repelem for high res
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        sqTmp = repelem(a,[round(sizeCycle(fq)*allDC(dc)) round(sizeCycle(fq)*(1-allDC(dc)))]);
%         figure;plot(sqTmp)
        % then copy it multiple times to get to the 1 sec window
        squareWave(:,fq,dc) = repmat(sqTmp,1,nT/sizeCycle(fq)*repCycle);
        
        idx1F1 =  (nT/sizeCycle(fq)*repCycle);
        
        idxAll = (idx1F1:idx1F1:(repCycle*nT/2) )+1;
        n = 1:length(idxAll);
        A = 2;
        d = allDC(dc);
        squareAnalyticFFT(idxAll,fq,dc) = abs(2 * A./(n*pi) .* sin ( n*pi*d));
        
        
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

maxFreq = 50;
figure; hold on
harm1 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(1)));
% plot(allDC,squeeze(squareFFT(harm1,1,1:end).^2),'b','LineWidth',2)
plot(allDC,sum(squeeze(squareFFT(harm1,1,1:end).^2)),'b','LineWidth',2)
f1=(sum(squeeze(squareFFT([idx1F1:idx1F1*2:maxFreq]+1,fq,:).^2)));
harm2 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(2)));
plot(allDC,squeeze(sum(squareFFT(harm2,2,1:end).^2)),'r','LineWidth',2)
harm3 = determineFilterIndices( 'nf1low49', freqs, find(freqs==allfq(3)));
plot(allDC,squeeze(sum(squareFFT(harm3,3,1:end).^2)),'g','LineWidth',2)
xlim([0 1])
legend('10Hz','5Hz','2Hz')
xlabel('duty cycle')
ylabel('sum2 fq amp < 50Hz')
title('normalisation')
saveas(gcf,['figures' filesep 'normalisationHighResCorrected'],'png')

figure; hold on
plot(allDC,sqrt(sum(squeeze(squareFFT(harm1,1,1:end).^2))),'b','LineWidth',2)
plot(allDC,sqrt(squeeze(sum(squareFFT(harm2,2,1:end).^2))),'r','LineWidth',2)
plot(allDC,sqrt(squeeze(sum(squareFFT(harm3,3,1:end).^2))),'g','LineWidth',2)
xlim([0 1])
legend('10Hz','5Hz','2Hz')
xlabel('duty cycle')
ylabel('sqrt(sum2 fq amp < 50Hz)')
title('normalisation')


%% Look at total  power (^2)

for iFreq  = 1:3,
fq = iFreq;
idx1F1 =  (nT/sizeCycle(fq)*repCycle);

% maxFreq = 100; 
%maxFreq = idx1F1*4+1;  %First 4 harmonics. 
%Use all freq
maxFreq = size(squareFFT,1); 

f1=(sum(squeeze(squareFFT([idx1F1:idx1F1*2:maxFreq]+1,fq,:).^2))); % odd harmonics
f2 =( sum(squeeze(squareFFT([idx1F1*2:idx1F1*2:maxFreq]+1,fq,:).^2))); % even harmonics
f1Analytic =( sum(squeeze(squareAnalyticFFT([idx1F1:idx1F1*2:maxFreq]+1,fq,:).^2)));
f2Analytic = ( sum(squeeze(squareAnalyticFFT([idx1F1*2:idx1F1*2:maxFreq]+1,fq,:).^2)));
 

figure

set(gcf,'DefaultLineLineWidth',2)
plot(allDC,f1)
hold on
plot(allDC,f1Analytic,'--')


plot(allDC,f2)
plot(allDC,f2Analytic,'--')

plot(allDC,(f1+f2)/2)
plot(allDC,(f1Analytic+f2Analytic)/2,'--')

title (['Power (sum of amp^2)'  num2str(allfq(fq)) ' Hz'])


%  Total Amplitude.  



f1=sum((squeeze(squareFFT([idx1F1:idx1F1*2:maxFreq]+1,fq,:))))
f2 = ( sum((squeeze(squareFFT([idx1F1*2:idx1F1*2:maxFreq]+1,fq,:)))))

f1Analytic = sum((squeeze(squareAnalyticFFT([idx1F1:idx1F1*2:maxFreq]+1,fq,:))))
f2Analytic = ( sum((squeeze(squareAnalyticFFT([idx1F1*2:idx1F1*2:maxFreq]+1,fq,:)))))
 

figure
title ('Power')
set(gcf,'DefaultLineLineWidth',2)
plot(allDC,f1)
hold on
plot(allDC,f1Analytic,'--')
plot(allDC,f2)

plot(allDC,f2Analytic,'--')

plot(allDC,(f1+f2)/2)
plot(allDC,(f1Analytic+f2Analytic)/2,'--')

title (['Sum of Amplitudes (sum of amp) ' num2str(allfq(fq)) ' Hz'])

end