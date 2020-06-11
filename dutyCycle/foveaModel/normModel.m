
load('minModel.mat');
load('16sbjDC.mat')

% remove 0 in the matrix
[nbElec, tot] = size(minGpModelRef.amp);
tmpAmp = nonzeros(getfield(minGpModelRef,'amp'));
tmpAmpReshape = reshape(tmpAmp,nbElec-1,length(tmpAmp)/(nbElec-1));
tmpCos = nonzeros(getfield(minGpModelRef,'cos'));
tmpCosReshape = reshape(tmpCos,nbElec-1,length(tmpCos)/(nbElec-1));

minModel.amp = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.amp(2:end,:) = tmpAmpReshape;
minModel.cos = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.cos(2:end,:) = tmpCosReshape;
% amplitude = abs(cos), no need to get it from matrix
% also there is no phase so sin = zeros
minModel.sin = zeros(size(minModel.cos));

% do a normalisation for each ROI -> unit norming
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
    if cond>15
        betaCos(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1*2-1), model);
        betaSin(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1*2-1), model);
        betaCosNoIPS(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1*2-1), model(:,[1:6 8]));
        betaSinNoIPS(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1*2-1), model(:,[1:6 8]));
    else
    betaCos(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1), model);
    betaSin(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1), model);
    betaCosNoIPS(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1), model(:,[1:6 8]));
    betaSinNoIPS(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1), model(:,[1:6 8]));
    end
end

betaAmp = sqrt(betaCos.^2 + betaSin.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmp(:,cond))
    ylim([0 25])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end

betaAmpNoIPS = sqrt(betaCosNoIPS.^2 + betaSinNoIPS.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmpNoIPS(:,cond))
    ylim([0 25])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','LOC'})
end
saveas(gcf,'betaNormalised','png')

figure;set(gcf, 'Position', [0 0 1500 600])
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    bar(betaAmpNoIPS(:,cond))
    ylim([0 25])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaNormalisedMotion','png')



