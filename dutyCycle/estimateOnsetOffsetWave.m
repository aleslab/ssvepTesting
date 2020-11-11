%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate onset, offset waveform and residual from static
% extrapole to motion signal
% then from the extracted waveforms (betas), fit/create new response for
% static and motion at different duty-cycles (designMatrix*betas) and give
% the results for 1f1 and per area (fit the area/eccentricity model)
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
%%
%%%% all freq together (stacked up)
maxnT = size(avData(15).wave,2);
timeLine = avData(15).time;
dutyCyclesPercent = repmat(dutyCyclesPercent,1,3);

for aa = 1:15 % 1:5 = 10Hz, 6:10=5Hz, 11:15=2.5Hz
    listnT(aa) = size(avData(aa).wave,2);
    dutyCyclesSamples(aa) = size(avData(aa).wave,2)*dutyCyclesPercent(aa);
    nbRep = maxnT/size(avData(aa).wave,2);
    data(:,:,aa) = repmat(avData(aa).wave,1,nbRep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% regression
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
% 
% for iDuty = 1:nDuty
% 
%     thisStack = (1:nT)+nT*(iDuty-1);%defines the stack location
%     designMatrix(thisStack,:)=  [eye(nT) circshift(eye(nT),dutyCyclesSamples(iDuty))];
% %     figure;imagesc(designMatrix(thisStack,:))
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
%%% Now let's do a quick plot to see how well the 2 waveforms fit the data.

%Plot the raw data, fit, and residual
figure('Position', [0 0 1000 500]);
subplot(2,3,1:2); hold on;
iCh = 23;
plot(stackedData(iCh,:)); %
plot(fitData(iCh,:));
plot(residual(iCh,:))
legend('Data','Fit','Residual')
title('allStatic')

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

saveas(gcf,['figures' filesep 'extractOnOffallStatic'],'png')


%%%% plot for bbsrc grant
iCh = 23;
figure; 
subplot(2,3,1);plot(squeeze(betaWeights(iCh,1,:))','g'); title('Onset response');ylim([-3 3])
xlabel('time (ms)')
ylabel('EEG amplitude')
subplot(2,3,2);plot(squeeze(betaWeights(iCh,2,:))','g'); title('Offset response');ylim([-3 3])
subplot(2,3,4); hold on;
plot(stackedData(iCh,192*11+1:192*12)); 
plot(fitData(iCh,192*11+1:192*12)); 
plot(residual(iCh,192*11+1:192*12)); title('25% duty-cycle')
xlabel('time (ms)')
ylabel('EEG amplitude')
subplot(2,3,5); hold on;
plot(stackedData(iCh,192*12+1:192*13)); 
plot(fitData(iCh,192*12+1:192*13));
plot(residual(iCh,192*12+1:192*13)); title('50% duty-cycle')
subplot(2,3,6); hold on;
plot(stackedData(iCh,192*13+1:192*14)); 
plot(fitData(iCh,192*13+1:192*14)); 
plot(residual(iCh,192*13+1:192*14)); title('75% duty-cycle')
ylim([-3 3])
legend('response','prediction','residuals','location','best')



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



%% use betas to fit the motion condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fit and residuals for the motion conditions based on the static on/off
%%% responses for left and right responses separately

dcMotion = [.125 .25 .5 .75 .875 .5 .5];
% prepare the data for left and right motion separately
for aa = 1:length(dcMotion)
    listnTmotion(aa) = size(avData(aa+15).wave,2)/2; 
    dcSamplesMotion(aa) = size(avData(aa+15).wave,2)/2*dcMotion(aa);
    nbRep = maxnT/(size(avData(aa+15).wave,2)/2);
    dataMotion(:,:,aa) = repmat(avData(aa+15).wave(:,1:end/2),1,nbRep);
    dataMotionLeft(:,:,aa) = repmat(avData(aa+15).wave(:,end/2+1:end),1,nbRep);
end
% create design matrix with hte motion condition and stacked motion data
for iCond = 1:length(dcMotion)
    thisStack = (1:nT)+nT*(iCond-1);
    matrixOn = repmat(eye(listnTmotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
    matrixOff = repmat(circshift(eye(listnTmotion(iCond)),dcSamplesMotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
    designMatrixMotion(thisStack,:) = [matrixOn matrixOff];
    for iCh = 1:nCh
        stackedDataMotion(iCh,thisStack) =  dataMotion(iCh,:,iCond);
        stackedDataMotionL(iCh,thisStack) =  dataMotionLeft(iCh,:,iCond);
    end
end
% use the beta found for the static conditions to create the fit for the
% motion data + residuals
for iCh = 1:nCh
    thisBetaMotion = [squeeze(betaWeights(iCh,1,:)); squeeze(betaWeights(iCh,2,:))]; 
    fitDataMotion(iCh,:) = (designMatrixMotion*thisBetaMotion)';
    residualMotion(iCh,:) = stackedDataMotion(iCh,:) - fitDataMotion(iCh,:);
    residualMotionL(iCh,:) = stackedDataMotionL(iCh,:) - fitDataMotion(iCh,:);
end

% %Plot the raw data, fit, and residual
% % the betas are from the static conditions, this is mostly for checks
% figure('Position', [0 0 1000 500]);
% subplot(2,3,1:2); hold on;
% iCh = 23;
% plot(stackedDataMotion(iCh,:)); %
% plot(fitDataMotion(iCh,:));
% plot(residualMotion(iCh,:))
% legend('Data','Fit','Residual')
% subplot(2,3,3)
% plot(squeeze(betaWeights(iCh,:,:))'); hold on;
% %1st and 2nd waveforms might correspond to "onset" and "offset" depending
% %on how the stimulus timing is coded.
% legend('1st waveform','2nd waveform')
% % plot topographies
% cfg.layout = 'biosemi128.lay';
% cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
% colormap('hot');
% subplot(2,3,4); plotTopo(rms(betaWeights(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
% subplot(2,3,5); plotTopo(rms(betaWeights(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
% subplot(2,3,6); plotTopo(rms(residual(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);
% saveas(gcf,['figures' filesep 'extractOnOffpredMotionRight'],'png')

% plot the fit for left and right motion with the residuals to compare
figure('Position', [0 0 1000 1000]);
subplot(3,3,1:3); hold on;
iCh = 23;
% plot data and fit
plot(stackedDataMotion(iCh,:)); %
plot(fitDataMotion(iCh,:));
plot(stackedDataMotionL(iCh,:)); %
legend('DataLeft','Fit','DataRight')
% plot residuals
subplot(3,3,4:6); hold on;
plot(residualMotion(iCh,:))
plot(residualMotionL(iCh,:))
legend('residualLeft','residualRight')
% plot topographies
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
colormap('hot');
subplot(3,3,7); plotTopo(rms(residualMotion(:,:),2),cfg.layout); title('residual right'); colorbar;caxis([0 0.6]);
subplot(3,3,8);plotTopo(rms(residualMotionL(:,:),2),cfg.layout); title('residual left'); colorbar;caxis([0 0.6]);
subplot(3,3,9);plotTopo(rms(residual(:,:),2),cfg.layout); title('residual static (left)'); colorbar;caxis([0 0.6]);
% ATTENTION: the stimulus starts on the left position but because I remove
% 3 (or multiple of 3) cycles, the waveform for the moving stimulus is
% first for a right stimulus then for a left stimulus
saveas(gcf,['figures' filesep 'extractOnOffpredMotion'],'png')

figure;colormap('hot');
subplot(1,2,1)
plotTopo(rms(residualMotionL,2) - rms(residual,2),cfg.layout); colorbar;caxis([0 0.2]);
title('residual diff left')
subplot(1,2,2)
plotTopo(rms(residualMotion,2) - rms(residual,2),cfg.layout); colorbar;caxis([0 0.2]);
title('residual diff right (static = left)')
saveas(gcf,['figures' filesep 'extractOnOffresidualsDiff'],'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% try using the entire motion window - should give same results
dcSamplesMotion = dcSamplesMotion*2;
nTmax = nT*2;
for aa = 1:length(dcMotion)
    nbRep = nTmax/(size(avData(aa+15).wave,2));
    dataMotionAll(:,:,aa) = repmat(avData(aa+15).wave,1,nbRep);
end
designMatrixMotAll=zeros(nTmax*length(dcMotion),nTmax*4);
for iCond = 1:length(dcMotion)
    thisStack = (1:nTmax)+nTmax*(iCond-1);
    matrixOnLeft = zeros(nTmax,nTmax); matrixOnRight = zeros(nTmax,nTmax);
    matrixOffLeft = zeros(nTmax,nTmax); matrixOffRight = zeros(nTmax,nTmax);
    
    matrixOnRight(1:nT,1:nT) = repmat(eye(listnTmotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
    matrixOnLeft(nT+1:end,nT+1:end) = repmat(eye(listnTmotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
    matrixOffRight(1:nT,1:nT) = repmat(circshift(eye(listnTmotion(iCond)),dcSamplesMotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
    matrixOffLeft(nT+1:end,nT+1:end) = repmat(circshift(eye(listnTmotion(iCond)),dcSamplesMotion(iCond)),nT/listnTmotion(iCond),nT/listnTmotion(iCond));
   
    designMatrixMotAll(thisStack,:) = [matrixOnRight matrixOnLeft matrixOffRight matrixOffLeft];
    for iCh = 1:nCh
        stackedDataAllMot(iCh,thisStack) =  dataMotionAll(iCh,:,iCond);
    end
end
for iCh = 1:nCh
    thisBetaMotion = zeros(nTmax*4,1);
    % add on response
    thisBetaMotion(1:nT) = squeeze(betaWeights(iCh,1,:));
    thisBetaMotion(nTmax+nT+1:nTmax*2) = squeeze(betaWeights(iCh,1,:));
    % add off response
    thisBetaMotion(nTmax*2+1:nTmax*2+nT) = squeeze(betaWeights(iCh,2,:));
    thisBetaMotion(nTmax*4-nT+1:end) = squeeze(betaWeights(iCh,2,:));
    fitDataMotAll(iCh,:) = (designMatrixMotAll*thisBetaMotion);
    residualMotAll(iCh,:) = stackedDataAllMot(iCh,:) - fitDataMotAll(iCh,:);
end
% plot
figure; hold on;
iCh = 23;
plot(stackedDataAllMot(iCh,:)); %
plot(fitDataMotAll(iCh,:));
plot(residualMotAll(iCh,:)); %
legend('Data','Fit','Residual')

figure; 
subplot(2,2,1);imagesc([matrixOnRight matrixOnLeft matrixOffRight matrixOffLeft])
subplot(2,2,2);hold on;plot(thisBetaMotion)
title('beta')
subplot(2,2,3);plot([squeeze(betaWeights(iCh,1,:)); squeeze(betaWeights(iCh,2,:))])
subplot(2,2,4);hold on;
plot(stackedDataAllMot(iCh,thisStack));
plot(fitDataMotAll(iCh,thisStack));
legend('data','fit')
saveas(gcf,['figures' filesep 'extractOnOffwhatswrongAgain'],'png')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute FFT for static and motion which will be use (cos and sin) to fit the area model
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
waveFormM = zeros(nCh,nT,7,2);
dataFFTM = zeros(nCh,length(freqs),7,2);
clear dataMotion
for cond = 1:7
    for ch=1:nCh
        clear dataFourierM dataFourierResM;
        dataMotion = fitDataMotion(ch,nT*(cond-1)+1:nT*cond); 
        waveFormM(ch,:,cond,1) = dataMotion;
        dataFourierM = fft(dataMotion);
        % normalise the Fourier (double freq in the signal and div by total
        % freq)
        dataFFTM(ch,:,cond,1) = 2*dataFourierM(1:nT/2+1)/nT;

        %         bar(freq,ampFFTnorm(ch,:,cond))
        dataResM = residualMotionL(ch,nT*(cond-1)+1:nT*cond);
        waveFormM(ch,:,cond,2) = dataResM;
        dataFourierResM = fft(dataResM);
        dataFFTM(ch,:,cond,2) = 2*dataFourierResM(1:nT/2+1)/nT;
    end
end

%% plot waveforms and freq for all DC all freq for both the fitdata and the residuals
% %%% figures
% allDC = [0.125 0.25 0.5 0.75 0.875];
% fqHz = [10 5 2];
% maxFreq = 20;
% tt = {'signal','residual'};
% 
% %%% Oz waveform + amplitude spectrum
% iChan = 23;
% for fqPres = 1:3
%     for type=1:2 % plot fit data =1 or residual =2
%         figure('Position', [0 0 1000 500]); hold on;
%         for dc=1:5
%             subplot(2,5,dc)
%             plot(timeLine,waveForm(iChan,:,dc+(fqPres-1)*5,type));
%             xlabel('Time (ms)')
%             ylabel('Amplitude');
%             title(['Oz ' num2str(fqHz(fqPres)) 'Hz ' tt{type} ' ' num2str(allDC(dc)*100)])
%             
%             subplot(2,5,dc+5)
%             bar(freqs(2:maxFreq),abs(dataFFT(iChan,2:maxFreq,dc+(fqPres-1)*5,type)));
%             xlabel('Frequency (Hertz)')
%             ylabel('Amplitude');
%         end
%         saveas(gcf,['figures' filesep 'onOffWave' num2str(fqHz(fqPres)) '-' num2str(type)],'png')
%     end
% end
% %%% topographies (rms of waveform)
% for type=1:2 % plot fit data =1 or residual =2
%     figure('Position', [0 0 1500 1000]); hold on;
%     for fqPres = 1:3
%         for dc=1:5
%             subplot(3,5,dc+(fqPres-1)*5)
%             plotTopo(rms(waveForm(:,:,dc+(fqPres-1)*5,type),2),cfg.layout);
%             colorbar;
%             title([num2str(fqHz(fqPres)) 'Hz ' tt{type} ' ' num2str(allDC(dc)*100)])
%         end
%     end
%     saveas(gcf,['figures' filesep 'onOffRMS-' num2str(type)],'png')
% end
%%
%%% the fitted data is constructed from the betaweights. The weights 
%%% represent the regression (common signal for onset/offset) across
%%% conditions. So the fitted data is created from a common onset/offset
%%% and simulated for the different duty-cycles using the design matrix. 
position = [11 7 3 9 15 13 8];
for type=1:2 % plot fit data =1 or residual =2
    figure('Position', [0 0 1500 1000]); hold on;
    for cond=1:size(waveFormM,3)
        subplot(3,5,position(cond))
        plotTopo(rms(waveFormM(:,:,cond,type),2),cfg.layout);
        colorbar;
    end
    saveas(gcf,['figures' filesep 'onOffRMSmotion-' num2str(type)],'png')
end
% difference static motion
figure('Position', [0 0 1500 1000]); hold on;
for cond=1:size(waveFormM,3)
    subplot(3,5,position(cond))
    plotTopo(rms(waveForm(:,:,position(cond),2),2) - rms(waveFormM(:,:,cond,2),2),cfg.layout);
    colorbar;
end
saveas(gcf,['figures' filesep 'onOffRMSdiffStaticMotion' ],'png')
%%
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the regression for 1f1 only
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
index1f1m = [2 3 5 3 2 2 3];
% motion
for type = 1:2
    for cond = 1:size(dataFFTM,3)
        betaCosM(:,cond,type) = regress (real(dataFFTM(:,index1f1m(cond),cond,type)) , model(:,[1:6 8])); % don't include IPS
        betaSinM(:,cond,type) = regress (-imag(dataFFTM(:,index1f1m(cond),cond,type)) , model(:,[1:6 8])); 
    end
end

betaAmp = sqrt(betaCos.^2 + betaSin.^2);
betaAmpM = sqrt(betaCosM.^2 + betaSinM.^2);

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
    ylabel('fitData amplitude 1f1')
end
saveas(gcf,['figures' filesep 'onOffbetaFitData1f1'],'png')


for type=1:2
    figure;set(gcf, 'Position', [0 0 1200 800])
    for area = 1:7
        subplot(2,4,area); hold on;
        plot(dcVal,betaAmp(area,1:5,type),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
        plot(dcVal,betaAmp(area,6:10,type),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
        plot(dcVal,betaAmp(area,11:15,type),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
        plot(dcVal(3),betaAmpM(area,3,type),['^:'  col{1}],'MarkerSize',15,'LineWidth',2)
        plot(dcVal([2 3 4]),betaAmpM(area,[2 7 4],type),['^:'  col{2}],'MarkerSize',15,'LineWidth',2)
        plot(dcVal([1 3 5]),betaAmpM(area,[1 6 5],type),['^:'  col{3}],'MarkerSize',15,'LineWidth',2)
        title(areaNames(area))
        xticks([0:12.5:100]);
        legend('10','5','2.5','Location','Best')
        xlabel('Duty Cycle')
        ylabel('fitData amplitude 1f1')
    end
    if type == 1
        saveas(gcf,['figures' filesep 'onOffbetaSignal'],'png')
    elseif type ==2
        saveas(gcf,['figures' filesep 'onOffbetaResiduals'],'png')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do the regression for all harmonics
% determine the indexes of all the harmonics < 50 Hz
% ATTENTION, need to use freqs which is a made-up list of harmonics
clear fitData fitDataM
divFq = repmat([8 16 32],5,1);divFq = divFq(:);
for cond=1:15
    tgtFq = 85/divFq(cond):85/divFq(cond):50;
    harmIdx{cond} = arrayfun(@(tgtFq) find(freqs == tgtFq),tgtFq);
end
divFqM = [32 16 8 16 32 32 16];
for cond=1:7
    tgtFq = 85/divFqM(cond):85/divFqM(cond):50;
    harmIdxM{cond} = arrayfun(@(tgtFq) find(freqs == tgtFq),tgtFq); 
end
% for static
for cond = 1:15
    for hh = 1 : length(harmIdx{cond})
        fitData(cond).betaCos(:,hh)  = regress (real(dataFFT(:,harmIdx{cond}(hh),cond,1)), model(:,[1:6 8]));
        fitData(cond).betaSin(:,hh)  = regress (-imag(dataFFT(:,harmIdx{cond}(hh),cond,1)), model(:,[1:6 8]));
        fitData(cond).betaCosRes(:,hh)  = regress (real(dataFFT(:,harmIdx{cond}(hh),cond,2)), model(:,[1:6 8]));
        fitData(cond).betaSinRes(:,hh)  = regress (-imag(dataFFT(:,harmIdx{cond}(hh),cond,2)), model(:,[1:6 8]));
    end
end
% motion
for cond = 1:7
    for hh = 1 : length(harmIdxM{cond})
        fitDataM(cond).betaCos(:,hh)  = regress (real(dataFFTM(:,harmIdxM{cond}(hh),cond,1)), model(:,[1:6 8]));
        fitDataM(cond).betaSin(:,hh)  = regress (-imag(dataFFTM(:,harmIdxM{cond}(hh),cond,1)), model(:,[1:6 8]));
        fitDataM(cond).betaCosRes(:,hh)  = regress (real(dataFFTM(:,harmIdxM{cond}(hh),cond,2)), model(:,[1:6 8]));
        fitDataM(cond).betaSinRes(:,hh)  = regress (-imag(dataFFTM(:,harmIdxM{cond}(hh),cond,2)), model(:,[1:6 8]));
    end
end

% compute the amplitude (for each of the harmonics)
% normalise with square function - attention indices for the harmonics are
% different for the square function and the computed dataFFT
for cond=1:length(avData)
   fitData(cond).harmIdx = determineFilterIndices( 'nf1low49', avData(cond).freq, avData(cond).i1f1);
end
for cond=1:7
   fitDataM(cond).harmIdx = determineFilterIndices( 'nf1low49', avData(cond+15).freq, avData(cond+15).i1f1);
end
load('squareFFT.mat') % load the POWER for a square function
for cond = 1:15
    % sum of power for the different areas
    fitData(cond).power = sum(fitData(cond).betaCos(:,:).^2 + fitData(cond).betaSin(:,:).^2,2);
    fitData(cond).powerRes = sum(fitData(cond).betaCosRes(:,:).^2 + fitData(cond).betaSinRes(:,:).^2,2);
    % amplitude
    fitData(cond).amp = sqrt(fitData(cond).power);
    fitData(cond).ampRes = sqrt(fitData(cond).powerRes);
    % normalised amplitudes
    sumSq(cond) = sum(sqAllFFT(avData(cond).harmIdx,cond));
    fitData(cond).ampNorm = sqrt(fitData(cond).power/ sumSq(cond));
    fitData(cond).ampNormRes = sqrt(fitData(cond).powerRes/ sumSq(cond));
end
for cond = 1:7
    fitDataM(cond).power = sum(fitDataM(cond).betaCos(:,:).^2 + fitDataM(cond).betaSin(:,:).^2,2);
    fitDataM(cond).powerRes = sum(fitDataM(cond).betaCosRes(:,:).^2 + fitDataM(cond).betaSinRes(:,:).^2,2);
    % amplitude
    fitDataM(cond).amp = sqrt(fitDataM(cond).power);
    fitDataM(cond).ampRes = sqrt(fitDataM(cond).powerRes);
    % normalised amplitudes
    sumSqM(cond) = sum(sqAllFFT(avData(cond+15).harmIdx,cond+15));
    fitDataM(cond).ampNorm = sqrt(fitDataM(cond).power/ sumSqM(cond));
    fitDataM(cond).ampNormRes = sqrt(fitDataM(cond).powerRes/ sumSqM(cond)); 
end



% easier format for plotting
plotSig = [fitData.ampNorm];
plotRes = [fitData.ampNormRes];
plotSigM = [fitDataM.ampNorm];
plotResM = [fitDataM.ampNormRes];

% plotSig = [fitData.amp];
% plotRes = [fitData.ampRes];
% plotSigM = [fitDataM.amp];
% plotResM = [fitDataM.ampRes];

% plot 7 regions * 22 conditions
col={'b','r','g'};
areaNames = {'V1','V2','V3','V3A','V4','MT','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(dcVal,plotSig(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,plotSig(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,plotSig(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,plotSigM(area,3),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],plotSigM(area,[2 7 4]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],plotSigM(area,[1 6 5]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(areaNames(area))
    xticks([0:12.5:100]);
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('fitData normAmp harmonics')
end
saveas(gcf,['figures' filesep 'onOffFitDataAllFqAmpNorm'],'png')

figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(dcVal,plotRes(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,plotRes(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,plotRes(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,plotResM(area,3),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],plotResM(area,[2 7 4]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],plotResM(area,[1 6 5]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(areaNames(area))
    xticks([0:12.5:100]);
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('fitData normAmp harmonics')
end
saveas(gcf,['figures' filesep 'onOffFitDataAllFqAmpNormRes'],'png')
%%

% % plot residual depending on ISI
% % this is because we thought that the residuals will be important to
% % explain percept but this is NOT the case!
% timeOn = timeLine(end) * (dcVal/100);
% timeOff = timeLine(end) - timeOn;
% 
% figure;set(gcf, 'Position', [0 0 1200 800])
% for area = 1:7
%     subplot(2,4,area); hold on;
%     plot(timeOff/4,betaAmp(area,1:5,2),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
%     plot(timeOff/2,betaAmp(area,6:10,2),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
%     plot(timeOff,betaAmp(area,11:15,2),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
%     title(areaNames(area))
%     xlim([0 400])
%     legend('10r','5r','2.5r','Location','Best')
%     xlabel('ISI')
%     ylabel('betaNormECC amplitude')
% end
% saveas(gcf,['figures' filesep 'onOffbeta-ISIresiduals'],'png')
% 
% figure;set(gcf, 'Position', [0 0 1200 800])
% for area = 1:7
%     subplot(2,4,area); hold on;
%     plot(timeOn/4,betaAmp(area,1:5,2),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
%     plot(timeOn/2,betaAmp(area,6:10,2),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
%     plot(timeOn,betaAmp(area,11:15,2),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
%     title(areaNames(area))
%     xlim([0 400])
%     legend('10r','5r','2.5r','Location','Best')
%     xlabel('On Time')
%     ylabel('betaNormECC amplitude')
% end


%% correlation signal and residuals with ratings ONLY FOR 1F1 
% ATTENTION: remember that the same beta is used to fit the static and the
% moving ratings!!!
load('fullRatings9.mat')
static = mean(tabStatic,3); 
flickRating = static(2:3,:)'; % do not include 10Hz
moving = mean(tabMot,3);
movRating = moving(2:3,:)';
movRating(movRating==0)=NaN;
tt = {'signal','residual'};
for type=1:2
    figure;set(gcf, 'Position', [0 0 1000 600])
    for area = 1:7
        subplot(2,4,area); hold on;
        scatter(betaAmp(area,6:15,type),(flickRating(:)), 80,'filled','MarkerEdgeColor','none');
        scatter(betaAmp(area,6:15,type),(movRating(:)), 80,'filled','MarkerEdgeColor','none');
        R = corrcoef(betaAmp(area,6:15,type),(flickRating(:)));
        RsqS = R(1,2).^2;
        Rm = corrcoef(betaAmp(area,[7 8 9 11 13 15],type),(movRating([2 3 4 6 8 10])));
        RsqM = Rm(1,2).^2;
        %     ylim([0 3]);
        xlabel(['beta from onOff' tt{type}])
        ylabel('motion rating')
        title([areaNames{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2)]) %only 2 digits
        lsline
    end
    saveas(gcf,['figures/onOffcorrelWithFit1f1' tt{type}],'png')
end


%% correlation signal and residuals with ratings for all freq static and motion
% ATTENTION: remember that the same beta is used to fit the static and the
% moving ratings so same x for static and moving
plotFit = [fitData(6:15).ampNorm]; % no 10 Hz
plotFitM = [fitDataM([1 2 4:7]).ampNorm];
figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(plotFit(area,:),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    scatter(plotFitM(area,:),moving([3 5 11 15 9 8]), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(plotFit(area,:),flickRating(:));
    RsqS = R(1,2).^2;
    Rm = corrcoef(plotFitM(area,:),(moving([3 5 11 15 9 8])));
    RsqM = Rm(1,2).^2;
    %     ylim([0 3]);
    xlabel('fit from static onOff norm Amp')
    ylabel('motion rating')
    title([areaNames{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2)]) %only 2 digits
    lsline
end
saveas(gcf,['figures/onOffcorrelWithFitAllFqNorm' ],'png')

plotFit = [fitData(6:15).ampNormRes]; % no 10 Hz
plotFitM = [fitDataM([1 2 4:7]).ampNormRes];
figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(plotFit(area,:),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    scatter(plotFitM(area,:),moving([3 5 11 15 9 8]), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(plotFit(area,:),flickRating(:));
    RsqS = R(1,2).^2;
    Rm = corrcoef(plotFitM(area,:),(moving([3 5 11 15 9 8])));
    RsqM = Rm(1,2).^2;
    %     ylim([0 3]);
    xlabel('fit from static onOff norm Amp')
    ylabel('motion rating')
    title([areaNames{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2)]) %only 2 digits
    lsline
end
saveas(gcf,['figures/onOffcorrelWithFitAllFqResNorm' ],'png')