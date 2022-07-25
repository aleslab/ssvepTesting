%%% short version sent to Justin

% load models
load('eccModel.mat');
% load data
load('16sbjDC.mat')



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
    if cond >15 % 1f1 becomes i1f1*2-1 for moving conditions
    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1*2-1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1*2-1), minModel.cos);
    else
    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1), minModel.cos);
    end
end

betaAmp = sqrt(betaCos.^2 + betaSin.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmp(:,cond))
    ylim([0 20])
    ylabel('betaAmp')
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaECC','png')



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

clear betaCos betaSin betaAmp
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

% with IPS
betaAmp = sqrt(betaCos.^2 + betaSin.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmp(:,cond))
    ylim([0 25])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end

% without IPS
% Cos
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaCosNoIPS(:,cond))
    ylim([-10 5])
    ylabel('betaCos')
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaNormCosECC','png')
% Sin
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaSin(:,cond))
    ylim([-10 15])
    ylabel('betaSinNoIPS')
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaNormSinECC','png')
% Amp
betaAmpNoIPS = sqrt(betaCosNoIPS.^2 + betaSinNoIPS.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmpNoIPS(:,cond))
    ylim([0 15])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','LOC'})
end
saveas(gcf,'betaNormECC','png')

figure;set(gcf, 'Position', [0 0 1500 600])
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    bar(betaAmpNoIPS(:,cond))
    ylim([0 25])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','LOC'})
end
saveas(gcf,'betaNormMotionECC','png')


% betaAmpNoIPS is 7 regions * 22 conditions
col={'b','r','g'};
tt = {'V1','V2','V3','V3A','V4','MT','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(dcVal,betaAmpNoIPS(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmpNoIPS(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,betaAmpNoIPS(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,betaAmpNoIPS(area,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],betaAmpNoIPS(area,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],betaAmpNoIPS(area,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(tt(area))
    xticks([0:12.5:100]);
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormECC amplitude')
end
saveas(gcf,'betaAmpPerArea','png')



