%%% normalisation of harmonics
% returns amplitudes at different fq for a square wave matching the 22
% stimulation conditions in the experiment (fq and duty cycle)
% need to square for power

clearvars;close all
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


%%%% create the response for all frequencies and duty-cycles
squareWave = zeros(length(t)*repCycle,length(allfq),length(allDC));

% number of time points per frequency cycle
for fq = 1:length(allfq)
    sizeCycle(fq) = length(0:sampling:allTimeWin(fq)-sampling);
end
a = [1 -1];
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        sqTmp = repelem(a,[round(sizeCycle(fq)*allDC(dc)) round(sizeCycle(fq)*(1-allDC(dc)))]);
        % then copy it multiple times to get to the 1 sec window
        squareWave(:,fq,dc) = repmat(sqTmp,1,nT/sizeCycle(fq)*repCycle);
    end
end


%%%% compute Fourier
squareFFT = zeros(length(freqs),length(allfq),length(allDC));
squareAnalyticFFT = squareFFT;
for fq = 1:length(allfq)
    for dc = 1 :length(allDC)
        clear tmpFFT
        tmpFFT = fft(squareWave(:,fq,dc));
        % normalise the Fourier
        squareFFT(:,fq,dc) = abs(2*tmpFFT(1:length(tmpFFT)/2+1)/length(tmpFFT));
        
        % compare with analytical FFT solution for square wave
        idx1F1 =  (nT/sizeCycle(fq)*repCycle);
        idxAll = (idx1F1:idx1F1:(repCycle*nT/2) )+1;
        n = 1:length(idxAll);
        A = 2;
        d = allDC(dc);
        squareAnalyticFFT(idxAll,fq,dc) = abs(2 * A./(n*pi) .* sin ( n*pi*d)); 
    end
end 


%%%% arrange for easy use in analysis
% need to permute so that when pooling fq & DC, it matches with
% the 22 conditions in the experiment
sqTmp = permute(squareFFT,[1 3 2]);
sqAtpm = permute(squareAnalyticFFT,[1 3 2]);
sqFFT = reshape(sqTmp,[size(sqTmp,1),size(sqTmp,2)*size(sqTmp,3)]);
sqAnalytic = reshape(sqAtpm,[size(sqAtpm,1),size(sqAtpm,2)*size(sqAtpm,3)]);

% add moving conditions (just keep the same normalisation) 
% 16:22 (fq DC): 2.5 12.5; 5 25; 10 50; 5 75; 2.5 87.5; 2.5 50; 5 50
eqCond = [11 7 3 9 15 13 8]; % corresponding static conditions for 16:22
for cond=16:22
    sqFFT(:,cond) =  sqFFT(:,eqCond(cond-15));
    sqAnalytic(:,cond) =  sqAnalytic(:,eqCond(cond-15));
end
save('squareFFTnew.mat','sqFFT','sqAnalytic')

% for fq = 1:length(allfq)
%     for dc = 1 :length(allDC)
%         idx1F1 =  (nT/sizeCycle(fq)*repCycle);
%         odd = sum(sqAnalytic(idx1F1+1:idx1F1*2:end,1));
%         even = sum(sqAnalytic(idx1F1*2+1:idx1F1*2:end,1));
%         (odd+even) / sum(sqAnalytic(idx1F1+1:idx1F1:end,1))
%     end
% end 

