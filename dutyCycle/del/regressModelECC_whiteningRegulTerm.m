% There is no hint of a U shape so noise covariance just makes things worse
% whatever regulterm is used and even with scaling factor

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

% determine the indexes of all the harmonics < 50 Hz
for cond=1:length(avData)
   harm(cond).idx = determineFilterIndices( 'nf1low49', avData(cond).freq, avData(cond).i1f1);
end

% plot the covariance for all noise, signal at 1f1, signal for all
% harmonics
noiseData = []; sigData = []; sigDataAll=[];
for cond = 1:length(avData)
    noiseData = [noiseData avData(cond).sin(:,harm(cond).idx-1) avData(cond).cos(:,harm(cond).idx-1)...
        avData(cond).sin(:,harm(cond).idx+1) avData(cond).cos(:,harm(cond).idx+1)];
end
for cond = 1:length(avData)
    sigData = [sigData avData(cond).sin(:,avData(cond).i1f1) avData(cond).cos(:,avData(cond).i1f1)...
        avData(cond).sin(:,avData(cond).i1f1) avData(cond).cos(:,avData(cond).i1f1)];
    sigDataAll = [sigDataAll avData(cond).sin(:,harm(cond).idx) avData(cond).cos(:,harm(cond).idx)...
        avData(cond).sin(:,harm(cond).idx) avData(cond).cos(:,harm(cond).idx)];
end
% compute covariance
noiseCov = cov(noiseData'); 
sigCov = cov(sigData'); sigCovAll = cov(sigDataAll'); 
figure; subplot(2,2,3);imagesc(noiseCov); colorbar; title('noise')
subplot(2,2,1);imagesc(sigCov); colorbar; title('1f1')
subplot(2,2,2);imagesc(sigCovAll); colorbar; title('sig harm <50Hz')
saveas(gcf,['figures' filesep 'covariancePlot'],'png')

% we need to find the right regularisation term. For this we will do a
% leave-one-out cross-validation: do the prewhitening on all electrodes
% except one, predict the signal of that electrode and compare with the
% recorded signal (compute least square error). Then sum all the sse and
% plot sum sse depending on the regulation parameter that was used. SSE
% should decrease with better fitting up to a point where it then increases
% due to underfitting (close to regul term 0 is overfitting: explains
% everything)
regulTerm = 1:2:26;
errorCov = zeros(length(regulTerm),1);
errorScale = zeros(length(regulTerm),1);
% errorId = zeros(length(regulTerm));
fitDiff = zeros(avData(1).nchan,length(avData));
for regul = 1:length(regulTerm)
    fprintf(['regulTerm ' num2str(regulTerm(regul)) '\n'])
    fitDiffCov = zeros(nbElec,length(avData));
    fitDiffScale = zeros(nbElec,length(avData));
%     fitDiffId = zeros(nbElec,length(avData));
    for rejElec = 1:nbElec
        chanSelect = 1:nbElec;
        chanSelect(chanSelect == rejElec) = [];
        
        % First collect all the "noise" 
        noiseData = [];
        for cond = 1:length(avData)
            noiseData = [noiseData avData(cond).sin(chanSelect,harm(cond).idx-1) avData(cond).cos(chanSelect,harm(cond).idx-1)...
                avData(cond).sin(chanSelect,harm(cond).idx+1) avData(cond).cos(chanSelect,harm(cond).idx+1)];
        end
        % compute covariance
        noiseCov = cov(noiseData');
        
        % noise all due to only a few components. Use regularisation 
        % parameter to increase the diagonal so that nothing is divided by 0
        noiseCov = noiseCov + eye(size(noiseCov))*max(diag(noiseCov))*regulTerm(regul);
        % figure; imagesc(noiseCov);
        
        % compare with a noiseCov based on identity matrix (ie each channel
        % only covaries with itself). Gives the same as regression without
        % covariance!
%         noiseId = eye(nbElec-1)*regulTerm(regul);
        
        % model fitting
        % lscov warning due to the reference channel at 0 which happens
        % with regulterm = 0
        betaCos = []; betaSin=[];betaCosCov = []; betaSinCov=[];
%         betaCosId = []; betaSinId=[];
        for cond = 1:length(avData)
            [betaCosCov(:,cond)]= lscov (minModel.cos(chanSelect,:),avData(cond).cos(chanSelect,avData(cond).i1f1), noiseCov,'orth');
            [betaSinCov(:,cond)]= lscov (minModel.cos(chanSelect,:),avData(cond).sin(chanSelect,avData(cond).i1f1), noiseCov,'orth');
%             [betaCosId(:,cond)]= lscov (minModel.cos(chanSelect,:),avData(cond).cos(chanSelect,avData(cond).i1f1), noiseId,'orth');
%             [betaSinId(:,cond)]= lscov (minModel.cos(chanSelect,:),avData(cond).sin(chanSelect,avData(cond).i1f1), noiseId,'orth');
            if regul == 1 % only needs to be done once (no change with regulterm)
                betaCos(:,cond)= regress (avData(cond).cos(chanSelect,avData(cond).i1f1), minModel.cos(chanSelect,:));
                betaSin(:,cond)= regress (avData(cond).sin(chanSelect,avData(cond).i1f1), minModel.cos(chanSelect,:));
            end
        end
        
        % predict the data from the betas = apply the betas to the model 
        % with all the electrodes
        dataFitCos = []; dataFitSin=[];dataFitCosCov = []; dataFitSinCov=[];
%         dataFitCosId = []; dataFitSinId=[];
        for cond = 1:length(avData)
%             dataFitCos(:,cond) = sum(betaCos(:,cond) .* minModel.cos');
%             dataFitSin(:,cond) = sum(betaSin(:,cond) .* minModel.cos');
            if regul == 1
                dataFitCos(:,cond) = minModel.cos * betaCos(:,cond);
                dataFitSin(:,cond) = minModel.cos * betaSin(:,cond);
            end
            dataFitCosCov(:,cond) = minModel.cos * betaCosCov(:,cond);
            dataFitSinCov(:,cond) = minModel.cos * betaSinCov(:,cond);
%             dataFitCosId(:,cond) = minModel.cos * betaCosId(:,cond);
%             dataFitSinId(:,cond) = minModel.cos * betaSinId(:,cond);
        end
%         test = sqrt(dataFitCos.^2 +  dataFitSin.^2);
%         figure;plotTopo(squeeze(test(:,8)),cfg.layout)
        
        % calculate the error for the left-out electrode
        % I guess should not do to all electrodes since they were put in at
        % the beginning? YEP
        for cond = 1:length(avData)
            if regul == 1
                fitDiff(rejElec,cond) = (avData(cond).cos(rejElec,avData(cond).i1f1) - dataFitCos(rejElec,cond))^2 +...
                    (avData(cond).sin(rejElec,avData(cond).i1f1) - dataFitSin(rejElec,cond))^2;
            end
            fitDiffCov(rejElec,cond) = (avData(cond).cos(rejElec,avData(cond).i1f1) - dataFitCosCov(rejElec,cond) )^2 +...
                (avData(cond).sin(rejElec,avData(cond).i1f1) - dataFitSinCov(rejElec,cond))^2;
%             %             fitDiffId(rejElec,cond) = (dataFitCosId(rejElec,cond) - avData(cond).cos(rejElec,avData(cond).i1f1))^2 +...
%             %                 (dataFitSinId(rejElec,cond) - avData(cond).sin(rejElec,avData(cond).i1f1))^2;
            
%             % look for a scale factor pb
%             scale = [avData(cond).cos(:,avData(cond).i1f1)' avData(cond).sin(:,avData(cond).i1f1)']...
%                 / [dataFitCosCov(:,cond)' dataFitSinCov(:,cond)'];
%             fitDiffScale(rejElec,cond) = (avData(cond).cos(rejElec,avData(cond).i1f1) - dataFitCosCov(rejElec,cond)*scale )^2 +...
%                 (avData(cond).sin(rejElec,avData(cond).i1f1) - dataFitSinCov(rejElec,cond)*scale)^2;
%             scaleDD(regul,rejElec,cond) = scale;
            
%             scaleCos = avData(cond).cos(:,avData(cond).i1f1)' / dataFitCosCov(:,cond)';
%             scaleSin = avData(cond).sin(:,avData(cond).i1f1)' / dataFitSinCov(:,cond)';
%             fitDiffScale(rejElec,cond) = (avData(cond).cos(rejElec,avData(cond).i1f1) - dataFitCosCov(rejElec,cond)*scaleCos )^2 +...
%                 (avData(cond).sin(rejElec,avData(cond).i1f1) - dataFitSinCov(rejElec,cond)*scaleSin)^2;
%             scaleCosDD(regul,rejElec,cond) = scaleCos;
%             scaleSinDD(regul,rejElec,cond) = scaleSin;
        end
%         % all elec
%         for cond = 1:length(avData)
%             if regul == 1
%                 fitDiffAll(rejElec,cond) = sum((dataFitCos(:,cond) - avData(cond).cos(:,avData(cond).i1f1)).^2 +...
%                     (dataFitSin(:,cond) - avData(cond).sin(:,avData(cond).i1f1)).^2);
%             end
%             fitDiffAllCov(rejElec,cond) = sum((dataFitCosCov(:,cond) - avData(cond).cos(:,avData(cond).i1f1)).^2 +...
%                 (dataFitSinCov(:,cond) - avData(cond).sin(:,avData(cond).i1f1)).^2);
% %             fitDiffId(rejElec,cond) = (dataFitCosId(rejElec,cond) - avData(cond).cos(rejElec,avData(cond).i1f1))^2 +...
% %                 (dataFitSinId(rejElec,cond) - avData(cond).sin(rejElec,avData(cond).i1f1))^2;
%         end
        
    end
    errorCov(regul) = sum(sum(fitDiffCov));
    errorScale(regul) = sum(sum(fitDiffScale));
    
%     errorAllCov(regul) = sum(sum(fitDiffAllCov));
%     errorId(regul) = sum(sum(fitDiffId));
end
errorReg = sum(sum(fitDiff)); % not taking care of noise covariance
% errorAllReg = sum(sum(fitDiffAll)); % not taking care of noise covariance

figure;hold on;
plot(regulTerm,errorCov)
line([regulTerm(1) regulTerm(end)],[errorReg errorReg],'Color','red')
% ylim([80 95])
% xlim([5 45])
saveas(gcf,['figures' filesep 'regulTerm'],'png')

% figure;hold on;
% plot(regulTerm,errorAllCov)
% line([regulTerm(1) regulTerm(end)],[errorAllReg errorAllReg],'Color','red')

figure;hold on;
plot(regulTerm,errorCov,'Color','b')
plot(regulTerm,errorScale,'Color','g')
line([regulTerm(1) regulTerm(end)],[errorReg errorReg],'Color','r')
legend('cov','cov+scale','no cov')
saveas(gcf,['figures' filesep 'regulTermScaled2'],'png')

figure;hold on;
for regul=1:9%length(regulTerm)
    subplot(3,3,regul)
    imagesc(squeeze(scaleDD(regul,:,:)))
    colorbar
    title(['Scale Regul' num2str(regulTerm(regul))])
end
saveas(gcf,['figures' filesep 'scalFact'],'png')
% presumably, if you add scale into the no cov for each condition it will 
% be at the green asymptote. There is no hint of a U shape so noise
% covariance just makes things worse

figure;hold on;
for regul=1:12%length(regulTerm)
    subplot(4,3,regul)
    imagesc(squeeze(scaleCosDD(regul,:,:)))
    colorbar
    title(['CosScale Regul' num2str(regulTerm(regul))])
end
saveas(gcf,['figures' filesep 'scalFactCos'],'png')

figure;hold on;
for regul=1:12%length(regulTerm)
    subplot(4,3,regul)
    imagesc(squeeze(scaleSinDD(regul,:,:)))
    colorbar
    title(['SinScale Regul' num2str(regulTerm(regul))])
end
saveas(gcf,['figures' filesep 'scalFactSin'],'png')
