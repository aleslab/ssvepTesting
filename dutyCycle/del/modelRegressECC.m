%%% long version (see also regressModelECC)
% regression
% 1/ look at the topographies for the different sources
% 2/ sanity check: when using fullModel, should find sources for a stim presented ventral + dorsal on the right visual field
% 3/ Dsin = ROI .* Betasin and Dcos = ROI .* Betacos
% BetaAmp = sqrt(Bsin2 + Bcos2)

addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
ft_defaults



% load models
load('eccModel.mat');
% load data
load('16sbjDC.mat')



% remove 0 in the matrix
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
%%%%%%%%%%% Look at the topographies for the different sources
%%%%%%%%%%% source = left ventral dorsal, ref = Cz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x y width height])
% locNames = {'V1','V2','V3','V3A','V4','MT','IPS','LOC'};
% figure('position', [200, 0, 1500, 800])
% colormap('hot')
% for src = 1:size(minModel.amp,2)
%     subplot(2,4,src)
%     plotTopo(squeeze(minModel.amp(:,src)),cfg.layout)
%     title(locNames{src})
% end

locNames = {'V1ecc','V2ecc','V3ecc','V3A','V4','MT','IPS','LOC'};
figure('position', [200, 0, 1500, 800])
colormap('parula')
for src = 1:size(minModel.amp,2)
    subplot(2,4,src)
    plotTopo(squeeze(minModel.cos(:,src)),cfg.layout)
    title(locNames{src})
end
saveas(gcf,'sourceLocCosECC','png')

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

% A column of 1's basically fits an average term and subtracts it out. 
% Basically it's transforming the data to average reference.  
% R^2 magnitude and f statistics are biased when an average term is present. 
% For most cases it makes sense to remove. Doesn?t matter for us at the moment.  
% We want a more strict fit. If we care we would do average reference.

for cond = 1:length(avData)
    if cond >15 % 1f1 becomes i1f1*2-1
    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1*2-1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1*2-1), minModel.cos);
    else
    [betaCos(:,cond),bintCos(:,:,cond),rCos(:,cond),rintCos(:,:,cond),statsCos(:,cond)]= regress (avData(cond).cos(:,avData(cond).i1f1), minModel.cos);
    [betaSin(:,cond),bintSin(:,:,cond),rSin(:,cond),rintSin(:,:,cond),statsSin(:,cond)]= regress (avData(cond).sin(:,avData(cond).i1f1), minModel.cos);
    end
end

figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaCos(:,cond))
    ylim([-10 10])
    ylabel('betaCos')
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
% saveas(gcf,'betaCosECC','png')

figure;set(gcf, 'Position', [0 0 1500 600])
for cond = 1:15
    subplot(3,5,cond)
    bar(betaSin(:,cond))
    ylim([-10 20])
    ylabel('betaSin')
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end
% saveas(gcf,'betaSinECC','png')

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

figure;set(gcf, 'Position', [0 0 1500 600])
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    bar(betaAmp(:,cond))
    ylim([0 20])
    set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','IPS','LOC'})
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% model without IPS
% %%% does not change the results for the other areas
% for cond = 1:length(avData)
%     betaCos2(:,cond)= regress (avData(cond).cos(:,avData(cond).i1f1), minModel.cos(:,[1:6 8]));
%     betaSin2(:,cond)= regress (avData(cond).sin(:,avData(cond).i1f1), minModel.cos(:,[1:6 8]));
% end
% 
% betaAmp2 = sqrt(betaCos2.^2 + betaSin2.^2);
% figure;set(gcf, 'Position', [0 0 1500 600])
% for cond = 1:15
%     subplot(3,5,cond)
%     bar(betaAmp2(:,cond))
%     ylim([-10 10])
%     set(gca, 'XTickLabel',{'V1','V2','V3','V3A','V4','MT','LOC'})
% end
% saveas(gcf,'betaAmp_noIPS','png')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create topographies from betas and model
clear totAmp;
for cond = 1: length(avData)
    totCos(:,cond) = minModel.cos * betaCos(:,cond);
    totSin(:,cond) = minModel.cos * betaSin(:,cond);
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
        caxis([0 2.5])
    end
    colorbar
end
saveas(gcf,'topoPredictedECC','png')


figure('position', [200, 0, 1500, 800])
colormap('hot')
position = [11 7 3 9 15 13 8];
for cond = 16:22
    subplot(3,5,position(cond-15))
    plotTopo(squeeze(totAmp(:,cond)),cfg.layout)
    if cond-15 == 3
        caxis([0 1.5]);
    else
        caxis([0 2.5]);
    end
    colorbar
end
saveas(gcf,'topoPredictedmotionECC','png')






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
%     ylim([0 15])
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormECC amplitude')
end
saveas(gcf,['figures' filesep 'betaAmpPerArea'],'png')







%%%% correlation with behavioural percept
load('fullRatings9.mat')
static = mean(tabStatic,3); moving = mean(tabMot,3);

flickRating = static(2:3,:)';
movRating = moving([3 5 11 15 9 8]);

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(betaAmpNoIPS(area,6:15),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(betaAmpNoIPS(area,6:15),flickRating(:));
    RsqS = R(1,2).^2;
    scatter(betaAmpNoIPS(area,[16:17 19:22]),movRating, 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(betaAmpNoIPS(area,[16:17 19:22]),movRating);
    RsqM = R(1,2).^2;
    R = corrcoef([betaAmpNoIPS(area,6:15) betaAmpNoIPS(area,[16:17 19:22])],[flickRating(:)' movRating]);
    RsqF = R(1,2).^2;
    ylim([0 3]);
    xlabel('betaNormECC')
    ylabel('motion rating')
    title([tt{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ' F=' num2str(RsqF,2)]) %only 2 digits
    lsline
end
% Rsquare reported for static, moving and both (F) conditions pooled
saveas(gcf,'figures/correlWithBeta','png')

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(betaAmpNoIPS(area,6:15),log(flickRating(:)), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(betaAmpNoIPS(area,6:15),log(flickRating(:)));
    RsqS = R(1,2).^2;
    scatter(betaAmpNoIPS(area,[16:17 19:22]),log(movRating), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(betaAmpNoIPS(area,[16:17 19:22]),log(movRating));
    RsqM = R(1,2).^2;
    R = corrcoef([betaAmpNoIPS(area,6:15) betaAmpNoIPS(area,[16:17 19:22])],[flickRating(:)' log(movRating)]);
    RsqF = R(1,2).^2;
%     ylim([0 3]);
    xlabel('betaNormECC')
    ylabel('log(motion rating)')
    title([tt{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ' F=' num2str(RsqF,2)]) %only 2 digits
    lsline
end
% Rsquare reported for static, moving and both (F) conditions pooled
saveas(gcf,'figures/correlWithBetaLog','png')





fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[1,0,0],...
               'Upper',[Inf,20,20],...
               'StartPoint',[1 1 1]);
           
g1 = fittype('-a*log(x)+b');
f1 = fit(betaAmpNoIPS(1,6:15)',flickRating(:),g1);
xx = linspace(1,11,50);
figure;plot(betaAmpNoIPS(1,6:15)',flickRating(:),'o',xx,f1(xx),'r-');

g = fittype('a*exp(-x)+b');
f0 = fit(betaAmpNoIPS(1,6:15)',flickRating(:),g);
xx = linspace(1,5,50);
figure;plot(betaAmpNoIPS(1,6:15)',flickRating(:),'o',xx,f0(xx),'r-');




x = betaAmpNoIPS(1,6:15)';
y = flickRating(:);
fcn = @(b,x) b(1)*exp(-x)+b(2) ;
[B,fval] = fminsearch(@(b) norm(y - fcn(b,x)), ones(2,1));



figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(betaAmpNoIPS(area,6:15),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    f0 = fit(betaAmpNoIPS(area,6:15)',flickRating(:),g);    
    xx = linspace(1,max(betaAmpNoIPS(area,6:15)),50);
    plot(xx,f0(xx),'r-');
    scatter(betaAmpNoIPS(area,[16:17 19:22]),movRating, 80,'filled','MarkerEdgeColor','none');
    f0 = fit(betaAmpNoIPS(area,[16:17 19:22])',movRating',g);    
    xx = linspace(1,max(betaAmpNoIPS(area,[16:17 19:22])),50);
    plot(xx,f0(xx),'r-');
    title([tt{area}]) %only 2 digits
end
saveas(gcf,'figures/funcExp','png')

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(betaAmpNoIPS(area,6:15),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    f1 = fit(betaAmpNoIPS(area,6:15)',flickRating(:),g1);    
    xx = linspace(1,max(betaAmpNoIPS(area,6:15)),50);
    plot(xx,f1(xx),'r-');
    scatter(betaAmpNoIPS(area,[16:17 19:22]),movRating, 80,'filled','MarkerEdgeColor','none');
    f1 = fit(betaAmpNoIPS(area,[16:17 19:22])',movRating',g1);    
    xx = linspace(1,max(betaAmpNoIPS(area,[16:17 19:22])),50);
    plot(xx,f1(xx),'r-');
    title([tt{area}]) %only 2 digits
end
saveas(gcf,'figures/funcLog','png')
