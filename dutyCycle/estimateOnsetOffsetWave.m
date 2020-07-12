%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate onset, offset waveform and residual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
ft_defaults
ft_defaultscfg.layout = 'biosemi128.lay';

%% get the data and choose what to put in the regression
% load('dataOnOffExtract.mat')
clearvars
load('16sbjDC.mat')
dutyCyclesPercent = [.125 .25 .5 .75 .875];


% %%%% each freq separately
% bb = 0;
% for aa = 11:15 % 1:5 = 10Hz, 6:10=5Hz, 11:15=2.5Hz
%     bb= bb+1;
%     listnT(bb) = size(avData(aa).wave,2);
%     data(:,:,bb) = avData(aa).wave;
% end
% if aa == 5
%     tt = '10Hz';
% elseif aa == 10
%     tt = '5Hz';
% elseif aa == 15
%     tt = '2Hz';
% end
% dutyCyclesSamples = size(avData(aa).wave,2)*dutyCyclesPercent;

%%%% all freq together (stacked up)
maxnT = size(avData(15).wave,2);
timeLine = avData(15).time;
bb = 0;
dutyCyclesPercent = repmat(dutyCyclesPercent,1,3);

for aa = 1:15 % 1:5 = 10Hz, 6:10=5Hz, 11:15=2.5Hz
    bb= bb+1;
    listnT(bb) = size(avData(aa).wave,2);
    dutyCyclesSamples(bb) = size(avData(aa).wave,2)*dutyCyclesPercent(aa);
    nbRep = maxnT/size(avData(aa).wave,2);
    data(:,:,bb) = repmat(avData(aa).wave,1,nbRep);
end
tt = 'all';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regression

nCh = size(data,1);%We will treat each electrode separately
nT  = size(data,2);
nDuty = size(data,3);

%First build "design matrix" that links duty cycles to timelines
%We will "stack" each duty cycle condition
%I'm going to do this explicityly instead of rely on tricky reshapes and
%matlab indexing, or by growing something every loop.  Just so it's clear
%and explicit how these are built.  indexing can get confusing

%For the data that results in nCh x nT*nDuty,
stackedData = zeros(nCh,nT*nDuty);
%So for the design matrix rows we have number of rows equal to nT * nDuty
%For columns we have nT*2 for the 2 different waveforms: onset and offset
designMatrix = zeros(nT*nDuty,nT*2);

%Now we are going to make a matrix to index into the "lags"
%The "Onset" one is simple (I'm assuming onset always starts at t=1
%T1=1, T2 = 2, ... Tlast = nT.  This is the same for all the duty cycles.
%And can be represented by the identity matrix which maps 1-1 and 2-2 and
%so forth.
%For the Offset it is shifted in time from the offset so for 50%:
%T1 = nT/2, T2 = nT/2 +1 ... Tnt/2 = 1,  Tlast = (nT/2) - 1
%It also happens at different times in the frame so for the different duty
%cycles, so we shift the start time by an amount equal to the duty cycle.
%This can be represented by taking the identity matrix and circular
%shifting the contents so the mapping is shifted.
%
%Note: I always make off by 1 errors.  So it's important to check that the
%amount that is shifted is correct.  and not +/- 1 off.  Haven't checked
%yet.
% for checking can use circshift(1:192,dutyCyclesSamples(iDuty))

% for iDuty = 1:nDuty
%
%     thisStack = (1:nT)+nT*(iDuty-1);%defines the stack location
%     designMatrix(thisStack,:)=  [eye(nT) circshift(eye(nT),dutyCyclesSamples(iDuty))];
%     for iCh = 1:nCh,
%         stackedData(iCh,thisStack) =  data(iCh,:,iDuty);
%     end
%
% end

%%%% the following is updated to accomodate different design matrix
%%%% depending on the frequency (but more confusing)
for iDuty = 1:nDuty
    thisStack = (1:nT)+nT*(iDuty-1);%defines the stack location
    matrixOn = repmat(eye(listnT(iDuty)),nT/listnT(iDuty),nT/listnT(iDuty));
    matrixOff = repmat(circshift(eye(listnT(iDuty)),dutyCyclesSamples(iDuty)),nT/listnT(iDuty),nT/listnT(iDuty));
    designMatrix(thisStack,:) = [matrixOn matrixOff];
    for iCh = 1:nCh
        stackedData(iCh,thisStack) =  data(iCh,:,iDuty);
    end
end


%Now lets work out the overlap waveforms by using a linear regression:
%stackedData(iCh) = designMatrix*betaWeights
%betaWeights = stackedData/designMatrix' Note the transpose.
%Now since the matrix is not square and is "singlular" we have to use the
%pseudo inverse: pinv().  Or some other regression method.

betaWeights = zeros(nCh,2,nT);

for iCh = 1:nCh
    %Do the regression to find the weightings.
    thisBeta = stackedData(iCh,:)*pinv(designMatrix');
    %Let's calculate the model fit:
    fitData(iCh,:) = (designMatrix*thisBeta')';
    %Let's also calculate the residual.  That's interesting to see how well
    %the linear summation model fits the data.
    %Note the transposes.  Always tricky to get this correct.
    residual(iCh,:) = stackedData(iCh,:) - fitData(iCh,:);
    
    %Next pull apart the onset (1) from the offset(2) to index more easily
    betaWeights(iCh,1,:) = thisBeta((1:nT));
    betaWeights(iCh,2,:) = thisBeta((1:nT)+nT);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% summary plot
%% Now let's do a quick plot to see how well the 2 waveforms fit the data.

%Plot the raw data, fit, and residual
figure('Position', [0 0 1000 500]);
subplot(2,3,1:2); hold on;
iCh = 23;
plot(stackedData(iCh,:)); %
plot(fitData(iCh,:));
plot(residual(iCh,:))
legend('Data','Fit','Residual')
title(tt)

subplot(2,3,3)
plot(squeeze(betaWeights(iCh,:,:))'); hold on;
%1st and 2nd waveforms might correspond to "onset" and "offset" depending
%on how the stimulus timing is coded.
legend('1st waveform','2nd waveform')

% plot topographies
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
colormap('hot');
subplot(2,3,4); plotTopo(rms(betaWeights(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,5); plotTopo(rms(betaWeights(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,6); plotTopo(rms(residual(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);

saveas(gcf,['figures' filesep 'extractOnOff' tt],'png')




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % do not include 1st condition (10Hz 12.5DC)
% for iCh = 1:nCh
%     thisBetaAlt = stackedData(iCh,nT+1:end)*pinv(designMatrix(nT+1:end,:)');
%     fitDataAlt(iCh,:) = (designMatrix(nT+1:end,:)*thisBetaAlt')';
%     residualAlt(iCh,:) = stackedData(iCh,nT+1:end) - fitDataAlt(iCh,:);
%     betaWeightsAlt(iCh,1,:) = thisBetaAlt((1:nT));
%     betaWeightsAlt(iCh,2,:) = thisBetaAlt((1:nT)+nT);
% end
% 
% %Plot the raw data, fit, and residual
% figure('Position', [0 0 1000 500]);
% subplot(2,3,1:2); hold on;
% iCh = 23;
% plot(stackedData(iCh,nT+1:end)); %
% plot(fitDataAlt(iCh,:));
% plot(residualAlt(iCh,:))
% legend('Data','Fit','Residual')
% title(tt)
% 
% subplot(2,3,3)
% plot(squeeze(betaWeightsAlt(iCh,:,:))'); hold on;
% %1st and 2nd waveforms might correspond to "onset" and "offset" depending
% %on how the stimulus timing is coded.
% legend('1st waveform','2nd waveform')
% 
% % plot topographies
% cfg.layout = 'biosemi128.lay';
% cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
% colormap('hot');
% subplot(2,3,4); plotTopo(rms(betaWeightsAlt(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
% subplot(2,3,5); plotTopo(rms(betaWeightsAlt(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
% subplot(2,3,6); plotTopo(rms(residualAlt(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);
% 
% saveas(gcf,['figures' filesep 'extractOnOff' tt],'png')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot waveforms and freq for all DC all freq for both the fitdata and the residuals
nbCond = length(fitData)/nT;
fsample = 510;
nfr = floor(nT/2)+1; %Number of unique frequency in fourier transform. The +1 is for the DC component.
freqs = (fsample/2)*linspace(0,1,nfr);% freq values
% first compute amplitude for the frequencies
waveForm = zeros(nCh,nT,nbCond,2);
dataFFT = zeros(nCh,length(freqs),nbCond,2);
for cond = 1:nbCond
    for ch=1:nCh
        clear dataFourier dataFourierRes;
        data = fitData(ch,nT*(cond-1)+1:nT*cond);
        waveForm(ch,:,cond,1) = data;
        dataFourier = fft(data);
        % normalise the Fourier (double freq in the signal and div by total
        % freq)
        dataFFT(ch,:,cond,1) = 2*dataFourier(1:nT/2+1)/nT;

        %         bar(freq,ampFFTnorm(ch,:,cond))
        dataRes = residual(ch,nT*(cond-1)+1:nT*cond);
        waveForm(ch,:,cond,2) = dataRes;
        dataFourierRes = fft(dataRes);
        dataFFT(ch,:,cond,2) = 2*dataFourierRes(1:nT/2+1)/nT;
    end
end

%%% figures
allDC = [0.125 0.25 0.5 0.75 0.875];
fqHz = [10 5 2];
maxFreq = 20;
tt = {'signal','residual'};

%%% Oz waveform + amplitude spectrum
iChan = 23;
for fqPres = 1:3
    for type=1:2 % plot fit data =1 or residual =2
        figure('Position', [0 0 1000 500]); hold on;
        for dc=1:5
            subplot(2,5,dc)
            plot(timeLine,waveForm(iChan,:,dc+(fqPres-1)*5,type));
            xlabel('Time (ms)')
            ylabel('Amplitude');
            title(['Oz ' num2str(fqHz(fqPres)) 'Hz ' tt{type} ' ' num2str(allDC(dc)*100)])
            
            subplot(2,5,dc+5)
            bar(freqs(2:maxFreq),abs(dataFFT(iChan,2:maxFreq,dc+(fqPres-1)*5,type)));
            xlabel('Frequency (Hertz)')
            ylabel('Amplitude');
        end
        saveas(gcf,['figures' filesep 'onOffWave' num2str(fqHz(fqPres)) '-' num2str(type)],'png')
    end
end
%%% topographies (rms of waveform)
for type=1:2 % plot fit data =1 or residual =2
    figure('Position', [0 0 1500 1000]); hold on;
    for fqPres = 1:3
        for dc=1:5
            subplot(3,5,dc+(fqPres-1)*5)
            plotTopo(rms(waveForm(:,:,dc+(fqPres-1)*5,type),2),cfg.layout);
            colorbar;
            title([num2str(fqHz(fqPres)) 'Hz ' tt{type} ' ' num2str(allDC(dc)*100)])
        end
    end
    saveas(gcf,['figures' filesep 'onOffRMS-' num2str(type)],'png')
end

%%% the fitted data is constructed from the betaweights. The weights 
%%% represent the regression (common signal for onset/offset) across
%%% conditions. So the fitted data is created from a common onset/offset
%%% and simulated for the different duty-cycles using the design matrix. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% instead of picking and electrode, plot results by area (using the eccentricity model)
% this is also problematic: the regression with the model is done on 1f1,
% but what about the other frequencies???

% size(dataFFT) : 128 97 15 2 (chan, freq, cond, signal/residual)
% amp = abs (dataFFT)
% cos = real (dataFFT)
% sin = -imag(dataFFT)

% load model
load('eccModel.mat');

% remove 0 in the model matrix
[nbElec, tot] = size(eccModelRef.amp);
tmpAmp = nonzeros(getfield(eccModelRef,'amp'));
tmpAmpReshape = reshape(tmpAmp,nbElec-1,length(tmpAmp)/(nbElec-1));
tmpCos = nonzeros(getfield(eccModelRef,'cos'));
tmpCosReshape = reshape(tmpCos,nbElec-1,length(tmpCos)/(nbElec-1));

minModel.amp = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.amp(2:end,:) = tmpAmpReshape;
minModel.cos = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.cos(2:end,:) = tmpCosReshape;
% amplitude = abs(cos), no need to get it from matrix
% also there is no phase so sin = zeros
minModel.sin = zeros(size(minModel.cos));

% normalise the model
model = minModel.cos;
modelNorm = sqrt(sum(minModel.cos.^2,1));
model = bsxfun(@rdivide,model,modelNorm);

% do the regression
% find the harmonics index for the 3 different freq
index1f1 = repmat([find(freqs == 85/8) find(freqs == 85/16) find(freqs == 85/32)],5,1);
index1f1 = index1f1(:);
betaCos = zeros(7,15,2);betaSin = zeros(7,15,2); % 7 areas
for type = 1:2
    for cond = 1:15
        betaCos(:,cond,type) = regress (real(dataFFT(:,index1f1(cond),cond,type)) , model(:,[1:6 8])); % don't include IPS
        betaSin(:,cond,type) = regress (-imag(dataFFT(:,index1f1(cond),cond,type)) , model(:,[1:6 8])); 
    end
end
betaAmp = sqrt(betaCos.^2 + betaSin.^2);

% plot beta for signal and residual
col={'b','r','g'};
areaNames = {'V1','V2','V3','V3A','V4','MT','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(dcVal,betaAmp(area,1:5,1),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,6:10,1),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,11:15,1),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,1:5,2),['.--'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,6:10,2),['.--'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,11:15,2),['.--'  col{3}],'MarkerSize',40,'LineWidth',2)
    title(areaNames(area))
    xticks([0:12.5:100]);
    legend('10s','5s','2.5s','10r','5r','2.5r','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormECC amplitude')
end
saveas(gcf,['figures' filesep 'onOffbeta'],'png')


% plot residual depending on ISI

timeOn = timeLine(end) * (dcVal/100);
timeOff = timeLine(end) - timeOn;

figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(timeOff/4,betaAmp(area,1:5,2),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(timeOff/2,betaAmp(area,6:10,2),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(timeOff,betaAmp(area,11:15,2),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    title(areaNames(area))
    xlim([0 400])
    legend('10r','5r','2.5r','Location','Best')
    xlabel('ISI')
    ylabel('betaNormECC amplitude')
end
saveas(gcf,['figures' filesep 'onOffbeta-ISIresiduals'],'png')

figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(timeOn/4,betaAmp(area,1:5,2),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(timeOn/2,betaAmp(area,6:10,2),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(timeOn,betaAmp(area,11:15,2),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    title(areaNames(area))
    xlim([0 400])
    legend('10r','5r','2.5r','Location','Best')
    xlabel('On Time')
    ylabel('betaNormECC amplitude')
end

