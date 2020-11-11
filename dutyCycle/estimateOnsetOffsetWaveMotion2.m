%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate onset, offset waveform and residual
% double the time to check for motion
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

%%%% all freq together (stacked up)
maxnT = size(avData(15).wave,2)*2;
dutyCyclesPercent = repmat(dutyCyclesPercent,1,3);

for aa = 1:15 % 1:5 = 10Hz, 6:10=5Hz, 11:15=2.5Hz
    listnT(aa) = size(avData(aa).wave,2);
    dutyCyclesSamples(aa) = size(avData(aa).wave,2)*2*dutyCyclesPercent(aa);
    nbRep = maxnT/size(avData(aa).wave,2);
    data(:,:,aa) = repmat(avData(aa).wave,1,nbRep);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regression
nCh = size(data,1);%We will treat each electrode separately
nT  = size(data,2);
nCond = size(data,3);

%For the data that results in nCh x nT*nDuty,
stackedData = zeros(nCh,nT*nCond);
%So for the design matrix rows we have number of rows equal to nT * nDuty
%For columns we have nT*2 for the 2 different waveforms: onset and offset
designMatrix = zeros(nT*nCond,nT*2);

%%%% the following is updated to accomodate different design matrix
%%%% depending on the frequency (but more confusing)
for iDuty = 1:nCond
    thisStack = (1:nT)+nT*(iDuty-1);%defines the stack location
    matrixOn = repmat(eye(listnT(iDuty)),nT/listnT(iDuty),nT/listnT(iDuty));
    matrixOff = repmat(circshift(eye(listnT(iDuty)),dutyCyclesSamples(iDuty)),nT/listnT(iDuty),nT/listnT(iDuty));
    designMatrix(thisStack,:) = [matrixOn matrixOff];
    for iCh = 1:nCh
        stackedData(iCh,thisStack) =  data(iCh,:,iDuty);
    end
end


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

