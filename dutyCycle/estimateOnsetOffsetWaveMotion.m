%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate onset, offset waveform and residual
% only applied for the static signal
% motion signal would need 2xtimeWindow for a left and right response
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
dutyCyclesPercent = [.125 .25 .5 .75 .875 .5 .5];

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
maxnT = size(avData(16).wave,2);
timeLine = avData(16).time;

for aa = 1:7 
    listnT(aa) = size(avData(aa+15).wave,2);
    dutyCyclesSamples(aa) = size(avData(aa+15).wave,2)*dutyCyclesPercent(aa);
    nbRep = maxnT/size(avData(aa+15).wave,2);
    data(:,:,aa) = repmat(avData(aa+15).wave,1,nbRep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regression for 4 betas (on/off left/right)
nCh = size(data,1);%We will treat each electrode separately
nT  = size(data,2);
nCond = size(data,3);

%First build "design matrix" that links duty cycles to timelines
%We will "stack" each duty cycle condition
%I'm going to do this explicityly instead of rely on tricky reshapes and
%matlab indexing, or by growing something every loop.  Just so it's clear
%and explicit how these are built.  indexing can get confusing

%For the data that results in nCh x nT*nDuty,
stackedData = zeros(nCh,nT*nCond);
%So for the design matrix rows we have number of rows equal to nT * nDuty
%For columns we have nT*4 for the 4 different waveforms: onset and offset
%for left and right response
designMatrix = zeros(nT*nCond,nT*4);

%%%% the following is updated to accomodate different design matrix
%%%% depending on the frequency (but more confusing)
for iCond = 1:nCond
    thisStack = (1:nT)+nT*(iCond-1);%defines the stack location
    matrixOn = repmat(eye(listnT(iCond)),nT/listnT(iCond),nT/listnT(iCond));
    matrixOnR = circshift(matrixOn,nT/2);
    matrixOff = repmat(circshift(eye(listnT(iCond)),dutyCyclesSamples(iCond)),nT/listnT(iCond),nT/listnT(iCond));
    matrixOffR = circshift(matrixOff,nT/2);
    designMatrix(thisStack,:) = [matrixOn matrixOnR matrixOff matrixOffR];
    for iCh = 1:nCh
        stackedData(iCh,thisStack) =  data(iCh,:,iCond);
    end
end

%Now lets work out the overlap waveforms by using a linear regression:
%stackedData(iCh) = designMatrix*betaWeights
%betaWeights = stackedData/designMatrix' Note the transpose.
%Now since the matrix is not square and is "singlular" we have to use the
%pseudo inverse: pinv().  Or some other regression method.

betaWeights = zeros(nCh,4,nT);

for iCh = 1:nCh
    %Do the regression to find the weightings.
    thisBeta = stackedData(iCh,:)*pinv(designMatrix');
    %Let's calculate the model fit:
    fitData(iCh,:) = (designMatrix*thisBeta')';
    %Let's also calculate the residual.  That's interesting to see how well
    %the linear summation model fits the data.
    %Note the transposes.  Always tricky to get this correct.
    residual(iCh,:) = stackedData(iCh,:) - fitData(iCh,:);
    
    %Next pull apart the onset (1) from the offset to index more easily
    betaWeights(iCh,1,:) = thisBeta((1:nT));
    betaWeights(iCh,2,:) = thisBeta((1:nT)+nT);
    betaWeights(iCh,3,:) = thisBeta((1:nT)+2*nT);
    betaWeights(iCh,4,:) = thisBeta((1:nT)+3*nT);
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
title('allMotion')

subplot(2,3,3)
plot(squeeze(betaWeights(iCh,:,:))'); hold on;
%1st and 2nd waveforms might correspond to "onset" and "offset" depending
%on how the stimulus timing is coded.
legend('ONleft','ONright','OFFright','OFFleft')

% plot topographies
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
colormap('hot');
subplot(2,3,4); plotTopo(rms(betaWeights(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,5); plotTopo(rms(betaWeights(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,6); plotTopo(rms(residual(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);

saveas(gcf,['figures' filesep 'extractOnOffMotion'],'png')


%% regression for 2 betas (on/off regardless of position)
stackedData2 = zeros(nCh,nT*nCond);
designMatrix2 = zeros(nT*nCond,nT*2);
listnT2 = listnT/2;

for iCond = 1:nCond
    thisStack = (1:nT)+nT*(iCond-1);%defines the stack location
    matrixOn = repmat(eye(listnT2(iCond)),nT/listnT2(iCond),nT/listnT2(iCond));
    matrixOff = repmat(circshift(eye(listnT2(iCond)),dutyCyclesSamples(iCond)),nT/listnT2(iCond),nT/listnT2(iCond));
    figure;imagesc([matrixOn matrixOff]);
    designMatrix2(thisStack,:) = [matrixOn matrixOff];
    for iCh = 1:nCh
        stackedData2(iCh,thisStack) =  data(iCh,:,iCond);
    end
end

%Now lets work out the overlap waveforms by using a linear regression:
%stackedData(iCh) = designMatrix*betaWeights
%betaWeights = stackedData/designMatrix' Note the transpose.
%Now since the matrix is not square and is "singlular" we have to use the
%pseudo inverse: pinv().  Or some other regression method.

betaWeights2 = zeros(nCh,2,nT);

for iCh = 1:nCh
    %Do the regression to find the weightings.
    thisBeta2 = stackedData2(iCh,:)*pinv(designMatrix2');
    %Let's calculate the model fit:
    fitData2(iCh,:) = (designMatrix2*thisBeta2')';
    %Let's also calculate the residual.  That's interesting to see how well
    %the linear summation model fits the data.
    %Note the transposes.  Always tricky to get this correct.
    residual2(iCh,:) = stackedData2(iCh,:) - fitData2(iCh,:);
    
    %Next pull apart the onset (1) from the offset (2) to index more easily
    betaWeights2(iCh,1,:) = thisBeta2((1:nT));
    betaWeights2(iCh,2,:) = thisBeta2((1:nT)+nT);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% summary plot
%% Now let's do a quick plot to see how well the 2 waveforms fit the data.

%Plot the raw data, fit, and residual
figure('Position', [0 0 1000 500]);
subplot(2,3,1:2); hold on;
iCh = 23;
plot(stackedData2(iCh,:)); %
plot(fitData2(iCh,:));
plot(residual2(iCh,:))
legend('Data','Fit','Residual')
title('allMotion')

subplot(2,3,3)
plot(squeeze(betaWeights2(iCh,:,:))'); hold on;
%1st and 2nd waveforms might correspond to "onset" and "offset" depending
%on how the stimulus timing is coded.
legend('ON','OFF')

% plot topographies
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
colormap('hot');
subplot(2,3,4); plotTopo(rms(betaWeights2(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,5); plotTopo(rms(betaWeights2(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,6); plotTopo(rms(residual2(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);

saveas(gcf,['figures' filesep 'extractOnOffMotion2'],'png')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% instead of a beta of 384 points (repeated), do 192 points. This doubles
%% the number of data
doubleData(:,:,1:7) = data(:,1:end/2,:); 
doubleData(:,:,8:14) = data(:,end/2+1:end,:);
nTd = size(doubleData,2);
dcSamples = repmat(dutyCyclesSamples,[1 2]);
listnTd = repmat(listnT2,[1 2]);
for iCond = 1:size(doubleData,3)
    thisStack = (1:nTd)+nTd*(iCond-1);%defines the stack location
    matrixOn = repmat(eye(listnTd(iCond)),nTd/listnTd(iCond),nTd/listnTd(iCond));
    matrixOff = repmat(circshift(eye(listnTd(iCond)),dcSamples(iCond)),nTd/listnTd(iCond),nTd/listnTd(iCond));
    designMatrixDouble(thisStack,:) = [matrixOn matrixOff];
    for iCh = 1:nCh
        stackedDataDouble(iCh,thisStack) =  doubleData(iCh,:,iCond);
    end
end
betaWeightsDouble = zeros(nCh,2,nTd);

for iCh = 1:nCh
    %Do the regression to find the weightings.
    thisBetaD = stackedData2(iCh,:)*pinv(designMatrixDouble');
    %Let's calculate the model fit:
    fitDataD(iCh,:) = (designMatrixDouble*thisBetaD')';
    %Let's also calculate the residual.  That's interesting to see how well
    %the linear summation model fits the data.
    %Note the transposes.  Always tricky to get this correct.
    residualD(iCh,:) = stackedDataDouble(iCh,:) - fitDataD(iCh,:);
    
    %Next pull apart the onset (1) from the offset (2) to index more easily
    betaWeightsDouble(iCh,1,:) = thisBetaD((1:nTd));
    betaWeightsDouble(iCh,2,:) = thisBetaD((1:nTd)+nTd);
end

%Plot the raw data, fit, and residual
figure('Position', [0 0 1000 500]);
subplot(2,3,1:2); hold on;
iCh = 23;
plot(stackedDataDouble(iCh,:)); %
plot(fitDataD(iCh,:));
plot(residualD(iCh,:))
legend('Data','Fit','Residual')
title('allMotion')

subplot(2,3,3)
plot(squeeze(betaWeightsDouble(iCh,:,:))'); hold on;
%1st and 2nd waveforms might correspond to "onset" and "offset" depending
%on how the stimulus timing is coded.
legend('ON','OFF')

% plot topographies
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
colormap('hot');
subplot(2,3,4); plotTopo(rms(betaWeightsDouble(:,1,:),3),cfg.layout); title('onset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,5); plotTopo(rms(betaWeightsDouble(:,2,:),3),cfg.layout); title('offset (beta)'); colorbar; %caxis([0 1.2]);
subplot(2,3,6); plotTopo(rms(residualD(:,:),2),cfg.layout); title('residual'); colorbar; %caxis([0 1.2]);

saveas(gcf,['figures' filesep 'extractOnOffMotionDouble'],'png')



%%%%%%% 
maxnT = size(avData(15).wave,2);
timeLine = avData(15).time;
dutyCyclesPercent = repmat(dutyCyclesPercent,1,3);
for aa = 1:15 % 1:5 = 10Hz, 6:10=5Hz, 11:15=2.5Hz
    listnT(aa) = size(avData(aa).wave,2);
    dutyCyclesSamples(aa) = size(avData(aa).wave,2)*dutyCyclesPercent(aa);
    nbRep = maxnT/size(avData(aa).wave,2);
    dataStatic(:,:,aa) = repmat(avData(aa).wave,1,nbRep);
end
cond = [11 7 3 9 15 13 8];
for iCond = 1:length(cond)
    thisStack = (1:nTd)+nTd*(iCond-1);%defines the stack location
    for iCh = 1:nCh
        stackedDataStatic(iCh,thisStack) =  dataStatic(iCh,:,cond(iCond));
    end
end


figure; hold on;
plot(stackedDataDouble(23,1:end/2),'b'); plot(stackedDataDouble(23,end/2+1:end),'g');
plot(stackedDataStatic(23,:),'r')
legend('left','right','static');
title('Oz response (~200 timepoints for 7 cond)')
saveas(gcf,['figures' filesep 'diffResp'],'png')

