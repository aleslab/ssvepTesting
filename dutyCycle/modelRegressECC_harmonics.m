% regression with the area loc model fitting the harmonics separately
clearvars

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


% determine the indexes of all the harmonics < 50 Hz
for cond=1:length(avData)
    if cond < 16
        fitData(cond).harmIdx = determineFilterIndices( 'nf1low49', avData(cond).freq, avData(cond).i1f1);
    else
        fitData(cond).harmIdx = determineFilterIndices( 'nf1low49', avData(cond).freq, avData(cond).i1f1*2-1);
    end
end


% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
model = minModel.cos;
modelNorm = sqrt(sum(minModel.cos.^2,1));
model = bsxfun(@rdivide,model,modelNorm);

 
% do the regression without including IPS
% this is without normalisation (useful for plotting the predicted maps)
for cond = 1:length(avData)
    for hh = 1 : length(fitData(cond).harmIdx)
        fitData(cond).betaCosRaw(:,hh)  = regress (avData(cond).cos(:,fitData(cond).harmIdx(hh)), minModel.cos(:,[1:6 8]));
        fitData(cond).betaSinRaw(:,hh)  = regress (avData(cond).sin(:,fitData(cond).harmIdx(hh)), minModel.cos(:,[1:6 8]));
    end
end

% do the regression without including IPS with a normalised model
for cond = 1:length(avData)
    for hh = 1 : length(fitData(cond).harmIdx)
        fitData(cond).betaCos(:,hh)  = regress (avData(cond).cos(:,fitData(cond).harmIdx(hh)), model(:,[1:6 8]));
        fitData(cond).betaSin(:,hh)  = regress (avData(cond).sin(:,fitData(cond).harmIdx(hh)), model(:,[1:6 8]));
    end
end

% compute the amplitude (for each of the harmonics)
for cond = 1:length(avData)
    fitData(cond).amp = sqrt(fitData(cond).betaCos(:,:).^2 + fitData(cond).betaSin(:,:).^2);
    % sum the amplitudes for all the harmonics while keeping the areas
    % separate
    fitData(cond).sumAmp = sum(fitData(cond).amp,2);
    % normalise the amplitudes
    load('squareFFT.mat') % load the amplitudes for a square function
    sumSq = sum(sqAllFFT(fitData(cond).harmIdx,cond));
    fitData(cond).normAmp = fitData(cond).sumAmp / sumSq;    
end



% easier format for plotting
allNormAmp = [fitData.normAmp];

% plot 7 regions * 22 conditions
col={'b','r','g'};
tt = {'V1','V2','V3','V3A','V4','MT','LOC'};
dcVal = [12.5 25 50 75 87.5];
figure;set(gcf, 'Position', [0 0 1200 800])
for area = 1:7
    subplot(2,4,area); hold on;
    plot(dcVal,allNormAmp(area,1:5),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,allNormAmp(area,6:10),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,allNormAmp(area,11:15),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,allNormAmp(area,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],allNormAmp(area,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],allNormAmp(area,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    title(tt(area))
    xticks([0:12.5:100]);
    ylim([0 12])
    legend('10','5','2.5','Location','Best')
    xlabel('Duty Cycle')
    ylabel('betaNormAllFq')
end
saveas(gcf,['figures' filesep 'betaNormAllFqPerArea_v2'],'png')




% have a look at how much each harmonic is involved in the amplitude
aa = [1 6 11];
fqTitle = {'10Hz','5Hz','2Hz'};

for fq=1:3
    allFq = [fitData(aa(fq):aa(fq)+4).amp]; % crappy order: DC1, harm1:maxHarm then DC2, harm1:maxHarm
    allFq = reshape(allFq,[7 length(allFq)/5 5]);
    
    figure;set(gcf, 'Position', [0 0 1200 800])
    for area = 1:7
        subplot(2,4,area); hold on;
        for harm = 1:4
            plot(squeeze(allFq(area,harm,:)),'LineWidth',2)
        end
        title([tt{area} '-' fqTitle{fq}])
        legend('h1','h2','h3','h4','Location','Best')
    end
    saveas(gcf,['figures' filesep 'betaNormHarmFq' fqTitle{fq}],'png')
end




%%%% correlation with behavioural percept
load('fullRatings9.mat')
static = mean(tabStatic,3); moving = mean(tabMot,3);

flickRating = static(2:3,:)';
movRating = moving([3 5 11 15 9 8]);

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(allNormAmp(area,6:15),flickRating(:), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,6:15),flickRating(:));
    RsqS = R(1,2).^2;
    scatter(allNormAmp(area,[16:17 19:22]),movRating, 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,[16:17 19:22]),movRating);
    RsqM = R(1,2).^2;
    R = corrcoef([allNormAmp(area,6:15) allNormAmp(area,[16:17 19:22])],[flickRating(:)' movRating]);
    RsqF = R(1,2).^2;
    ylim([0 3]);
    xlabel('betaNormECC')
    ylabel('motion rating')
    title([tt{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ' F=' num2str(RsqF,2)]) %only 2 digits
%     lsline
end
% Rsquare reported for static, moving and both (F) conditions pooled
saveas(gcf,['figures' filesep 'correlWithBetaNorm'],'png')

figure;set(gcf, 'Position', [0 0 1000 600])
for area = 1:7
    subplot(2,4,area); hold on;
    scatter(allNormAmp(area,6:15),log(flickRating(:)), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,6:15),log(flickRating(:)));
    RsqS = R(1,2).^2;
    scatter(allNormAmp(area,[16:17 19:22]),log(movRating), 80,'filled','MarkerEdgeColor','none');
    R = corrcoef(allNormAmp(area,[16:17 19:22]),log(movRating));
    RsqM = R(1,2).^2;
    R = corrcoef([allNormAmp(area,6:15) allNormAmp(area,[16:17 19:22])],[flickRating(:)' log(movRating)]);
    RsqF = R(1,2).^2;
%     ylim([0 3]);
    xlabel('betaNormECC')
    ylabel('log(motion rating)')
    title([tt{area} ' S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ' F=' num2str(RsqF,2)]) %only 2 digits
%     lsline
end
% Rsquare reported for static, moving and both (F) conditions pooled
saveas(gcf,['figures' filesep 'correlWithBetaLogNorm'],'png')

