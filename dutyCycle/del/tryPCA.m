%%% PCA on sin cosin for all electrodes
% again, pca done on only the first harmonic for now

% load data
load('16sbjDC.mat')

% data should be in the form of 128 elect * conditions
dataIn = zeros(128,30);
for cond = 1:15
    dataIn(:,cond) = avData(cond).sin(:,avData(cond).i1f1);
    dataIn(:,cond+15) = avData(cond).cos(:,avData(cond).i1f1);
end

clear coeff score latent tsquared explained mu
[coeff,score,latent,tsquared,explained,mu] = pca(dataIn');

figure;set(gcf, 'Position', [0 0 1200 800])
for cpnt=1:6
    subplot(2,3,cpnt)
    plotTopo(coeff(:,cpnt),cfg.layout)
    colorbar
    title(num2str(explained(cpnt),'%.2f'))
end
saveas(gcf,['figures' filesep 'pca'],'png')

% i don't understand if using this one:
clear coeff score latent tsquared explained mu
[coeff,score,latent,tsquared,explained,mu] = pca(dataIn);
% get a 30x30 matrix of coefs
% explained: [78.52;13.29;5.43;1.29;0.72;0.23]

% pca using inverse of variances as variable weights = crap
clear coeff score latent tsquared explained mu
[coeff,score,latent,tsquared,explained,mu] = pca(dataIn(2:end,:)','VariableWeights','variance');

