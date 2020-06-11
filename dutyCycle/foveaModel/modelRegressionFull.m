% regression
% 1/ look at the topographies for the different sources
% 2/ sanity check: when using fullModel, should find sources for a stim presented ventral + dorsal on the right visual field
% 3/ Dsin = ROI .* Betasin and Dcos = ROI .* Betacos
% BetaAmp = sqrt(Bsin2 + Bcos2)

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults



% load models
load('fullModel.mat');
% load data
load('15sbjDC.mat')

% plot the actual (observed) topographies
figure('position', [200, 0, 1500, 800])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(squeeze(avData(cond).amp(:,avData(cond).i1f1)),cfg.layout)
%     colorbar;
    caxis([0 2.5])
end



% remove 0 in the matrix
[nbElec, tot] = size(fullGpModelRef.amp);
tmpAmp = nonzeros(getfield(fullGpModelRef,'amp'));
tmpAmpReshape = reshape(tmpAmp,nbElec-1,length(tmpAmp)/(nbElec-1));
tmpCos = nonzeros(getfield(fullGpModelRef,'cos'));
tmpCosReshape = reshape(tmpCos,nbElec-1,length(tmpCos)/(nbElec-1));

fullModel.amp = zeros(nbElec,length(tmpAmp)/(nbElec-1));
fullModel.amp(2:end,:) = tmpAmpReshape;
fullModel.cos = zeros(nbElec,length(tmpAmp)/(nbElec-1));
fullModel.cos(2:end,:) = tmpCosReshape;
fullModel.sin = zeros(size(fullModel.cos));

for cond = 1:15 % length(avData)
    [betaCosFull(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), fullModel.cos);
    [betaSinFull(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1), fullModel.cos);
end

betaAmpFull = sqrt(betaCosFull.^2 + betaSinFull.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
aa=0;
for cond = [1 2 3 6 7 8 11 12 13]
    aa = aa+1;
    subplot(3,3,aa)
    bar(betaAmpFull(:,cond))
%     ylim([0 15])
%     set(gca, 'XTickLabel',{'V1V-L','V1V-R','V1D-L','V1D-R',...
%     'V2V-L','V2V-R','V2D-L','V2D-R',...
%     'V3V-L','V3V-R','V3D-L','V3D-R','V3A-L','V3A-R',...
%     'V4-L','V4-R','MT-L','MT-R',...
%     'IPS-L','IPS-R','LOC-L','LOC-R'})
    set(gca, 'XTick', [1:2:22])
    set(gca, 'XTickLabel',{'V1V','V1D',...
    'V2V','V2D',...
    'V3V','V3D','V3A',...
    'V4','MT',...
    'IPS','LOC'})
end
saveas(gcf,'betaAmpFull','png')

for cond = 1: 15
    totCos = zeros(128,1);
    totAmp = zeros(128,1);
    for roi=1:size(betaAmpFull,1)
        totCos = totCos + betaAmpFull(roi,cond) * fullModel.cos(:,roi); % take cos or amp?? Seems not to matter
        totAmp = totAmp + betaAmpFull(roi,cond) * fullModel.cos(:,roi); % take cos or amp?? Seems not to matter
   end
    predictedTopoFullCos(:,cond) = totCos;
    predictedTopoFullAmp(:,cond) = totAmp;
end

figure('position', [200, 0, 1500, 800])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(squeeze(predictedTopoFullCos(:,cond)),cfg.layout)
%     colorbar;
%     caxis([0 10])
end
figure('position', [200, 0, 1500, 800])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(squeeze(predictedTopoFullAmp(:,cond)),cfg.layout)
%     colorbar;
    caxis([0 22])
end
saveas(gcf,'predictedTopoScaledFullAmp','png')

