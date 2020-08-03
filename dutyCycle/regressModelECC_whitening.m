% load models
load('eccModel.mat');
% load data
load('16sbjDC_v2.mat')

addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated

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

%Normalize column power
%Note! this uses implicit binary singleton expansion
minModel.cos = minModel.cos./sum(minModel.cos.^2);


%% %%
%%%% Prewhitening
%%%

% determine the indexes of all the harmonics < 50 Hz
for cond=1:length(avData)
   harm(cond).idx = determineFilterIndices( 'nf1low49', avData(cond).freq, avData(cond).i1f1);
end

% %First collect all the "noise"
% noiseData = [];
% sigData = [];
% for cond = 1:length(avData)
%     avData(cond).i1f1
%     nfMax = avData(cond).i1f1*4;
%     nf1= ((avData(cond).i1f1):(avData(cond).i1f1-1):nfMax)
%     noiseIdx = [nf1-1 nf1+1]
%     noiseData = [noiseData avData(cond).sin(:,noiseIdx) avData(cond).cos(:,noiseIdx)];
%     sigData = [sigData avData(cond).sin(:,nf1) avData(cond).cos(:,nf1)];
% end
% 
% %Given the noise values we calculate the noise only covariance matrix
% noiseCov = cov(noiseData');
% %"prewhitening" means removing the bias from the correlated noise.  
% %That means basically dividing the data by the covariance matrix. But since
% %EEG data is heavily correlated it
% W = cholcov(pinv(noiseCov));
% 
% sigCov = cov(sigData');
% figure;imagesc(sigCov);
% figure;imagesc(noiseCov);


chStart = 1; % 1=all channels, 2=do not include reference channel
% does not seem to matter at all so just keep it to avoid confusion
noiseData2 = [];
sigData2 = [];
for cond = 1:length(avData)
    sigData2 = [sigData2 avData(cond).sin(chStart:end,harm(cond).idx) avData(cond).cos(chStart:end,harm(cond).idx)];
    noiseData2 = [noiseData2 avData(cond).sin(chStart:end,harm(cond).idx-1) avData(cond).cos(chStart:end,harm(cond).idx-1)...
        avData(cond).sin(chStart:end,harm(cond).idx+1) avData(cond).cos(chStart:end,harm(cond).idx+1)];
end
noiseCov2 = cov(noiseData2');
sigCov2 = cov(sigData2');
figure;imagesc(sigCov2);colorbar;title('signal covariance');saveas(gcf,['figures' filesep 'covSig'],'png')
figure;imagesc(noiseCov2);colorbar;title('noise covariance');saveas(gcf,['figures' filesep 'covNoise'],'png')


% noise all due to only a few components
% regularisation parameter to increase the diagonal so that nothing is
% divided by 0
noiseCov = noiseCov2 + eye(size(noiseCov2))*max(diag(noiseCov2))*1;

[u,s,v]=svd(noiseCov);
figure;plot(diag(s))

%MVDR beamforming weight matrix.  
%mvdrW = pinv(noiseCov)*minModel.cos
%W= W';
%W=eye(128);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regression: Data = ROI * Betas. Matlab: Y = X*B, B = regress(Y,X)
% X is an n-by-p design matrix, with rows corresponding to observations and 
% columns to predictor variables
% Y is an n-by-1 vector of response observations.
% X should include a column of ones so that the model contains a constant term.
% Since I am looking at the response at the harmonics I think I should only
% use that frequency

% could use mldivide as well: give the same results
% M = rand(8,3); % model 8 electrodes 3 ROIS
% D = rand(8,1); % data 8 electrodes cosinus of 1 harmonic (1f1)
% beta = M\D;
% beta2 = regress(D,M);

for cond = 1:length(avData)
    noiseOff =0;
    [betaCos(:,cond)]= lscov (minModel.cos(chStart:end,:),avData(cond).cos(chStart:end,avData(cond).i1f1+noiseOff), noiseCov,'orth');
    [betaSin(:,cond)]= lscov (minModel.cos(chStart:end,:),avData(cond).sin(chStart:end,avData(cond).i1f1+noiseOff), noiseCov,'orth');
end

betaAmp = sqrt(betaCos.^2 + betaSin.^2);
% figure;set(gcf, 'Position', [0 0 1500 600])
% for cond = 1:15
%     subplot(3,5,cond)
%     bar(betaAmp(:,cond))
%     ylim([0 max(betaAmp(:))])
%     ylabel('betaAmp')
%     set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'},'fontsize',12)
% end
% 
% figure(101);
% for iRoi = 1:size(betaAmp,1)
%     
%     subplot(2,4,iRoi)
%     plot(betaAmp(iRoi,:));
%     hold on;
% end


col={'b','r','g'};
tt = {'V1','V2','V3','V3A','V4','MT','IPS','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:8
    subplot(2,4,area); hold on;
    plot(dcVal,betaAmp(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmp(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,betaAmp(area,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],betaAmp(area,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],betaAmp(area,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(tt(area))
    xticks([0:12.5:100]);
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormWhite amplitude')
end
saveas(gcf,['figures' filesep 'betaAmpPerAreaWhite'],'png')




%% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1

model = minModel.cos;

modelNorm = sqrt(sum(minModel.cos.^2,1));

model = bsxfun(@rdivide,model,modelNorm);

figure; imagesc(minModel.cos)
figure; imagesc(model)

for cond = 1:length(avData)
    noiseOff =0;
    betaCosNorm(:,cond)= lscov (model,avData(cond).cos(:,avData(cond).i1f1+noiseOff), noiseCov,'orth');
    betaSinNorm(:,cond)= lscov (model,avData(cond).sin(:,avData(cond).i1f1+noiseOff), noiseCov,'orth');
end

betaAmpNorm = sqrt(betaCosNorm.^2 + betaSinNorm.^2);
col={'b','r','g'};
tt = {'V1','V2','V3','V3A','V4','MT','IPS','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:8
    subplot(2,4,area); hold on;
    plot(dcVal,betaAmpNorm(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmpNorm(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmpNorm(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,betaAmpNorm(area,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],betaAmpNorm(area,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],betaAmpNorm(area,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(tt(area))
    xticks([0:12.5:100]);
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormWhite amplitude')
end
saveas(gcf,['figures' filesep 'betaAmpPerAreaNormWhite'],'png')



