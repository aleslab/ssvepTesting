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
load('minModel.mat');
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
    if cond <6
        caxis([0 1.5])
    else
        caxis([0 3])
    end
end
saveas(gcf,'topoDC','png')
% motion
figure('position', [200, 0, 1500, 800])
position = [11 7 3 9 15 13 8];
colormap('hot')
for cond=16:length(avData)
%     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avData(cond).amp(:,avData(cond).i1f1*2-1),cfg.layout);
    if cond-15 == 3
        caxis([0 1.5]);
    else
        caxis([0 3]);
    end
end
saveas(gcf,'topoDCmotion','png')


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Look at the topographies for the different sources
%%%%%%%%%%% source = left ventral dorsal, ref = Cz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height])
locNames = {'V1','V2','V3','V3A','V4','MT','IPS','LOC'};
figure('position', [200, 0, 1500, 800])
colormap('hot')
for src = 1:size(minModel.amp,2)
    subplot(2,4,src)
    plotTopo(squeeze(minModel.amp(:,src)),cfg.layout)
    title(locNames{src})
end
saveas(gcf,'sourceLoc','png')

locNames = {'V1','V2','V3','V3A','V4','MT','IPS','LOC'};
figure('position', [200, 0, 1500, 800])
colormap('parula')
for src = 1:size(minModel.amp,2)
    subplot(2,4,src)
    plotTopo(squeeze(minModel.cos(:,src)),cfg.layout)
    title(locNames{src})
end
saveas(gcf,'sourceLocCos','png')

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

% A column of 1?s basically fits an average term and subtracts it out. 
% Basically it?s transforming the data to average reference.  
% R^2 magnitude and f statistics are biased when an average term is present. 
% For most cases it makes sense to remove. Doesn?t matter for us at the moment.  
% We want a more strict fit. If we care we would do average reference.

for cond = 1:length(avData)
    if cond >15 % 1f1 becomes i1f1*2-1
    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1*2-1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1*2-1), minModel.cos);
        
    else
    % with intercept term
%     [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), [ones(length(minModel.cos),1) minModel.cos]);
%     [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), [ones(length(minModel.cos),1) minModel.cos]);

    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1), minModel.cos);
    % should give same results
    betaCos2(:,cond) = minModel.cos\avData(cond).cos(:,avData(cond).i1f1);
    betaSin2(:,cond) = minModel.cos\avData(cond).sin(:,avData(cond).i1f1);
    end
end

figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaCos(:,cond))
    ylim([-10 10])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaCos','png')
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaSin(:,cond))
    ylim([-10 10])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaSin','png')

betaAmp = sqrt(betaCos.^2 + betaSin.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmp(:,cond))
    ylim([0 15])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaAmp','png')

figure;set(gcf, 'Position', [0 0 1500 600])
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    bar(betaAmp(:,cond))
    ylim([0 15])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
saveas(gcf,'betaAmpMotion','png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% model without IPS
%%% does not change the results for the other areas
for cond = 1:length(avData)
    betaCos2(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1), minModel.cos(:,[1:6 8]));
    betaSin2(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1), minModel.cos(:,[1:6 8]));
end

betaAmp2 = sqrt(betaCos2.^2 + betaSin2.^2);
figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaAmp2(:,cond))
    ylim([-10 10])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','LOC'})
end
saveas(gcf,'betaAmp_noIPS','png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create topographies from betas and model
clear totAmp;
for cond = 1: length(avData)
    totCos(:,cond) = zeros(128,1);
    totSin(:,cond) = zeros(128,1);
    for roi=1:size(betaCos,1)
        totCos(:,cond) = totCos(:,cond) + betaCos(roi,cond) * minModel.cos(:,roi);
        totSin(:,cond) = totSin(:,cond) + betaSin(roi,cond) * minModel.cos(:,roi);
    end
    totAmp(:,cond) = sqrt(totCos(:,cond).^2 + totSin(:,cond).^2); 
end

figure('position', [200, 0, 1500, 800])
colormap('hot')
for cond = 1:15
    subplot(3,5,cond)
    plotTopo(squeeze(totAmp(:,cond)),cfg.layout)
    if cond <6
        caxis([0 1.5])
    else
        caxis([0 3])
    end
end
saveas(gcf,'topoPredicted','png')


figure('position', [200, 0, 1500, 800])
colormap('hot')
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    plotTopo(squeeze(totAmp(:,cond)),cfg.layout)
    if cond-15 == 3
        caxis([0 1.5]);
    else
        caxis([0 3]);
    end
end
saveas(gcf,'topoPredictedmotion','png')


