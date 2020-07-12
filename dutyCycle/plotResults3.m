
clearvars
% addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
% ft_defaults

% addpath /Users/Marlene/Documents/git/fieldtrip-aleslab-fork
% addpath /Users/Marlene/Documents/git/ssvepTesting/svndlCopy
% addpath /Users/Marlene/Documents/git/ssvepTesting/biosemiUpdated
% addpath /Users/Marlene/Documents/git/ssvepTesting/commonFunctions
% ft_defaults

addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
ft_defaults

% load individual data
% dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
% dataDir = 'C:\Users\Marlene\Documents\JUSTIN\data\dutyCycle\Axx\';
dataDir = '/Volumes/Amrutam/Marlene/JUSTIN/DutyCycle/data/Axx/';

listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
dcVal = [12.5 25 50 75 87.5];

keepSbj = [1:5 7:11 13:15 17:18 20]; 
% reject S7 and S13 S17, S20 same as S03
% S1=pilote so reject numS-1 (6 12 16 19)

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    dataSbj(:,ss) = Axx;
end

[avData, proj_Amp] = averageAxxWithStd(dataSbj);
save('16sbjDC','avData','proj_Amp','cfg')

    

col={'b','r','g'};

% pickElec = 23; % 23 Oz, 9, B6=38, 16=best SNR

% from topo max: A15 & A27
pickElec = [15 23 27];


%%% create a matrix for the harmonics depending on condition
harm = zeros(15,4);
for cond=1:15
    filtIdx = determineFilterIndices( 'nf1low49',avData(cond).freq, avData(cond).i1f1 );
    harm(cond,:) = filtIdx(1:4); % avData(1).freq(filtIdx)
end
% add the harmonics for motion stim
harm(16,:) = harm(11,:);
harm(21,:) = harm(11,:);
harm(20,:) = harm(11,:);
harm(18,:) = harm(1,:);
harm(17,:) = harm(6,:);
harm(19,:) = harm(6,:);
harm(22,:) = harm(6,:);


% 2.5 - 5 - 7.5 - 10
% 5 - 10 - 15 - 20
% 10 - 20 - 30 - 40

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Topography
% for each harmonic separately
% flicker
for hh = 1:4
    figure; hold on;
    colormap('hot')
    for cond=1:15
%         max(avData(cond).amp(:,harm(cond,hh)))
        subplot(3,5,cond); hold on;
        plotTopo(avData(cond).amp(:,harm(cond,hh)),cfg.layout);
        colorbar
        if hh == 1 && cond < 6
            caxis([0 1.5]);
        elseif hh == 1 && cond > 5
            caxis([0 2.5]);
        elseif hh == 2 && cond > 10
            caxis([0 2]);
        elseif hh == 2 && cond > 5 && cond < 11
            caxis([0 1]);
        elseif hh == 3 && cond > 5
            caxis([0 0.8]);
        elseif hh == 4 && cond > 10
            caxis([0 0.8]);
        else
            caxis([0 0.5]);
        end
        
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1500 1000])
    saveas(gcf,['figures' filesep 'topoFlickF' num2str(hh) '.jpg'])
end
% motion
position = [11 7 3 9 15 13 8];
for hh = 1:4
    figure; hold on;
    colormap('hot')
    for cond=16:length(avData)
%         max(avData(cond).amp(:,harm(cond,hh)))
        subplot(3,5,position(cond-15)); hold on;
        plotTopo(avData(cond).amp(:,harm(cond,hh)),cfg.layout);
        colorbar
        if hh > 1 && cond == 18
            caxis([0 0.5]);
        elseif hh == 1 && cond  == 18
            caxis([0 1.5]);
        elseif hh == 1 && cond  ~= 18
            caxis([0 2.5]);
        elseif hh == 3 && cond ~= 18
            caxis([0 0.8]);
        elseif hh == 2 && ismember(cond,[16 21 20])
            caxis([0 2]);
        elseif hh == 2 && ismember(cond,[17 22 19])
            caxis([0 1]);
        elseif hh == 4 && ismember(cond,[16 21 20])
            caxis([0 0.8]);
        else
            caxis([0 0.5]);
        end
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1500 1000])
    saveas(gcf,['figures' filesep 'topoMotionF' num2str(hh) '.jpg'])
end



% TOPO sum across harmonics 
for cond=1:15
    fqharm(cond).ind(:) = determineFilterIndices( 'nf1low49',avData(cond).freq, avData(cond).i1f1 );
end
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(sum(avData(cond).amp(:,fqharm(cond).ind),2),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 2]);
    elseif cond > 5 && cond < 11
        caxis([0 4.5]);
    elseif cond > 10
        caxis([0 6]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickSumFq.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    %     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(sum(avData(cond).amp(:,fqharm(position(cond-15)).ind),2),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 2]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 4.5]);
    elseif position(cond-15) > 10
        caxis([0 6]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionSumFq.jpg'])
% normalise the sum 
load squareFFT.mat
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(sum(avData(cond).amp(:,fqharm(cond).ind),2)/sum(sqAllFFT(fqharm(cond).ind,cond)),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1]);
    elseif cond > 5 && cond < 11
        caxis([0 2]);
    elseif cond > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickAllFqNorm.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
    %     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(sum(avData(cond).amp(:,fqharm(position(cond-15)).ind),2)/sum(sqAllFFT(fqharm(position(cond-15)).ind,cond)),cfg.layout);
    colorbar
    if position(cond-15) < 6
        caxis([0 1]);
    elseif position(cond-15) > 5 && position(cond-15) < 11
        caxis([0 2]);
    elseif position(cond-15) > 10
        caxis([0 2.5]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionAllFqNorm.jpg'])


%%%%%
normAmp = zeros(128,22);
for cond=1:15
    normAmp(:,cond) = sum(avData(cond).amp(:,fqharm(cond).ind),2)/sum(sqAllFFT(fqharm(cond).ind,cond));
end
for cond = 16:22
    normAmp(:,cond) = sum(avData(cond).amp(:,fqharm(position(cond-15)).ind),2)/sum(sqAllFFT(fqharm(position(cond-15)).ind,cond));
end

col={'b','r','g'};
chan = 23;
figure;hold on;
for freq=1:3
    plot(dcVal,normAmp(chan,(freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',40,'Linewidth',2);
end
% line([dcVal(1) dcVal(5)], [min(baseline) min(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
% line([dcVal(1) dcVal(5)], [mean(baseline) mean(baseline)],'Linewidth',2,'Color','k','LineStyle','-')
% line([dcVal(1) dcVal(5)], [max(baseline) max(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
plot(50,normAmp(chan,18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
plot([25 50 75],normAmp(chan,[17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
plot([12.5 50 87.5],normAmp(chan,[16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
legend('10Hz','5Hz','2.5Hz','Location','Best')
xlabel('Duty Cycle')
ylabel('norm SSVEP amplitude')
set(gca,'FontSize',15)
title('Oz')
saveas(gcf,['figures' filesep 'ampDCNormAllFq.jpg'])



%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SSVEP

% get F1 amplitude
for elec = 1:length(pickElec)
    for cond=1:15
        avAmp(cond)=avData(cond).amp(pickElec(elec),(avData(cond).i1f1));
        %    avProj(cond)=mean(proj_Amp(pickElec,avData(cond).i1f1,cond,:));
        semProj(cond) = std(proj_Amp(pickElec(elec),avData(cond).i1f1,cond,:)) / sqrt(length(avData));
    end
    for cond=16:22 % for motion take the 2f1
        avAmp(cond)=avData(cond).amp(pickElec(elec),(avData(cond).i1f1*2-1));
        semProj(cond) = std(proj_Amp(pickElec(elec),avData(cond).i1f1*2-1,cond,:)) / sqrt(length(avData));
    end
    
    % compute noise level
    for cond=1:15 
       baseline(cond) = (avData(cond).amp(pickElec(elec),avData(cond).i1f1-1)+avData(cond).amp(pickElec(elec),avData(cond).i1f1+1)) / 2;
    end
    for cond=16:22
       baseline(cond) = (avData(cond).amp(pickElec(elec),avData(cond).i1f1*2-1-1)+avData(cond).amp(pickElec(elec),avData(cond).i1f1*2-1+1)) / 2;
    end
    
    
    
    figure;hold on;
    for freq=1:3
        errorbar(dcVal,avAmp((freq-1)*5+1:freq*5),semProj((freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',40,'Linewidth',2);
    end
    line([dcVal(1) dcVal(5)], [min(baseline) min(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
    line([dcVal(1) dcVal(5)], [mean(baseline) mean(baseline)],'Linewidth',2,'Color','k','LineStyle','-')
    line([dcVal(1) dcVal(5)], [max(baseline) max(baseline)],'Linewidth',2,'Color','k','LineStyle','--')
    errorbar(50,avAmp(18),semProj(18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    errorbar([25 50 75],avAmp([17 22 19]),semProj([17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    errorbar([12.5 50 87.5],avAmp([16 21 20]),semProj([16 21 20]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    xlim([0 100]);
    xticks([0:12.5:100]);
%     xticklabels({'12.5','25','50','75','87.5'})
    ylim([0 3])
    legend('10','5','2.5','noise range','mean noise','Location','Best')
    xlabel('Duty Cycle')
    ylabel('SSVEP amplitude')
    set(gca,'FontSize',15)
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'ampDCloc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'ampDCloc' num2str(pickElec(elec)) '.eps'],'epsc')
    
    %%%%%% LINEAR REGRESSION
    %%%%%% this is on average data not each individual data points
    % fitlm(dcVal,avAmp(1:5))
    % fitlm(dcVal,avAmp(6:10))
    % fitlm(dcVal,avAmp(11:15))
    %
    % fitlm(dcVal(2:5),avAmp(2:5))
    
    
    % polyfit(dcVal,avAmp(6:10),1)
    % corrcoef(dcVal,avAmp(6:10))
    
    
    %%% OUTPUT:
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                    Estimate        SE        tStat     pValue
    %                    _________    _________    ______    _______
    %     (Intercept)      0.44001      0.30684     1.434    0.24703
    %     x1             0.0079936    0.0053312    1.4994    0.23073
    %
    % Number of observations: 5, Error degrees of freedom: 3
    % Root Mean Squared Error: 0.34
    % R-squared: 0.428,  Adjusted R-Squared 0.238
    % F-statistic vs. constant model: 2.25, p-value = 0.231
    % ans =
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                    Estimate        SE         tStat      pValue
    %                    _________    _________    _______    _________
    %     (Intercept)       2.7107      0.31251     8.6739    0.0032242
    %     x1             -0.028999    0.0054298    -5.3407     0.012835
    %
    % Number of observations: 5, Error degrees of freedom: 3
    % Root Mean Squared Error: 0.346
    % R-squared: 0.905,  Adjusted R-Squared 0.873
    % F-statistic vs. constant model: 28.5, p-value = 0.0128
    % ans =
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                    Estimate        SE        tStat       pValue
    %                    _________    _________    ______    __________
    %     (Intercept)       2.6201      0.19678    13.314    0.00091569
    %     x1             -0.025466    0.0034191    -7.448     0.0050102
    %
    % Number of observations: 5, Error degrees of freedom: 3
    % Root Mean Squared Error: 0.218
    % R-squared: 0.949,  Adjusted R-Squared 0.932
    % F-statistic vs. constant model: 55.5, p-value = 0.00501
    
    
    
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                     Estimate        SE         tStat      pValue
    %                    __________    _________    _______    ________
    %     (Intercept)       0.99393      0.15841     6.2744    0.024473
    %     x1             0.00026448    0.0024735    0.10693     0.92461
    %
    % Number of observations: 4, Error degrees of freedom: 2
    % Root Mean Squared Error: 0.119
    % R-squared: 0.00568,  Adjusted R-Squared -0.491
    % F-statistic vs. constant model: 0.0114, p-value = 0.925
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%% RATINGS
    load('fullRatings9.mat')
    static = mean(tabStatic,3); moving = mean(tabMot,3);
    for fq=1:3
        for dc=1:5
            statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(6);
            movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(6);
        end
    end
    
    if elec == 1 % don't plot for each electrode
        figure;hold on;
        for fq=1:3
            errorbar(dcVal,static(fq,:),statSEM(1,:),['.-'  col{fq}],'MarkerSize',40,'LineWidth',2)
        end
        errorbar(dcVal(3),moving(7),movSEM(7),['^'  col{1}],'MarkerSize',15,'Linewidth',2)
        errorbar([25 50 75],moving([5 8 11]),movSEM([5 8 11]),['^:'  col{2}],'MarkerSize',15,'Linewidth',2)
        errorbar([12.5 50 87.5],moving([3 9 15]),movSEM([3 9 15]),['^:'  col{3}],'MarkerSize',15,'Linewidth',2)
        xlim([0 100]);
        xticks([0:12.5:100]);
        ylim([0 3])
        % xticklabels({'12.5','25','50','75','87.5'})
        legend('10','5','2.5','Location','Best')
        xlabel('Duty Cycle')
        ylabel('motion ratings')
        set(gca,'FontSize',15)
        saveas(gcf,'ratingsDC.png')
        saveas(gcf,'ratingsDC.eps','epsc')
    end
    
    % for freq=1:3
    %     fitlm(dcVal,static(freq,1:5))
    % end
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%% CORREL
    y = [static(2,:) static(3,:)];
    yE = [statSEM(2,:) statSEM(3,:)];
    x = avAmp(6:15);
    xE = semProj(6:15);
    y2 = moving([3 5 11 15 9 8]); % do not include 10Hz
    yE2 = movSEM([3 5 11 15 9 8]);
    x2 = avAmp([16:17 19:22]);% do not include 10Hz
    xE2 = semProj([16:17 19:22]);
    
    % figure;hold on;
    % p = polyfit(x,y,1);
    % f = polyval(p,x);
    % plot(x,y,'.',x,f,'-')
    % p2 = polyfit(x2,y2,1);
    % f2 = polyval(p2,x2);
    % plot(x2,y2,'^',x2,f2,'-')
    % ylim([0 3])
    % legend('flickering','','moving')
    % ylabel('motion rating')
    % xlabel('SSVEP amplitude')
    % title(num2str(pickElec(elec)))
    % saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.png'])
    % saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.eps'],'epsc')
    

    figure;hold on;
    scatter(x, y,80,'filled','MarkerEdgeColor','none');
    scatter(x2, y2,80,'filled','MarkerEdgeColor','none');
    eb(1) = errorbar(x,y,xE, 'horizontal', 'LineStyle', 'none');
    eb(2) = errorbar(x,y,yE, 'vertical', 'LineStyle', 'none');
    set(eb, 'color', 'b', 'LineWidth', 1)
    eb2(1) = errorbar(x2,y2,xE2, 'horizontal', 'LineStyle', 'none');
    eb2(2) = errorbar(x2,y2,yE2, 'vertical', 'LineStyle', 'none');
    set(eb2, 'color', 'r', 'LineWidth', 1)
    lsline
    ylim([0 3])
    legend('flickering','moving')
    ylabel('motion rating')
    xlabel('SSVEP amplitude')
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.pdf'],'pdf')
  
    figure;hold on;
    scatter(log(x), log(y),80,'filled','MarkerEdgeColor','none');
    scatter(log(x2), log(y2),80,'filled','MarkerEdgeColor','none');
    lsline
    legend('flickering','moving')
    ylabel('motion rating (log scale)')
    xlabel('SSVEP amplitude (log scale)')
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLogLoc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLogLoc' num2str(pickElec(elec)) '.pdf'],'pdf')
      
    figure; hold on;
    scatter(x(1:5), y(1:5),80,'filled','MarkerFaceColor','b','MarkerEdgeColor','none');
    scatter(x(6:10), y(6:10),80,'d','filled','MarkerFaceColor','c','MarkerEdgeColor','none');
    scatter(x2(1:3), y2(1:3),80,'filled','MarkerFaceColor','r','MarkerEdgeColor','none');
    scatter(x2(4:6), y2(4:6),80,'d','filled','MarkerFaceColor','m','MarkerEdgeColor','none');
    eb(1) = errorbar(x(1:5),y(1:5),xE(1:5), 'horizontal', 'LineStyle', 'none');
    eb(2) = errorbar(x(1:5),y(1:5),yE(1:5), 'vertical', 'LineStyle', 'none');
    set(eb, 'color', 'b', 'LineWidth', 1)
    eb(1) = errorbar(x(6:10),y(6:10),xE(6:10), 'horizontal', 'LineStyle', 'none');
    eb(2) = errorbar(x(6:10),y(6:10),yE(6:10), 'vertical', 'LineStyle', 'none');
    set(eb, 'color', 'c', 'LineWidth', 1)
    eb2(1) = errorbar(x2(1:3),y2(1:3),xE2(1:3), 'horizontal', 'LineStyle', 'none');
    eb2(2) = errorbar(x2(1:3),y2(1:3),yE2(1:3), 'vertical', 'LineStyle', 'none');
    set(eb2, 'color', 'r', 'LineWidth', 1)
    eb2(1) = errorbar(x2(4:6),y2(4:6),xE2(4:6), 'horizontal', 'LineStyle', 'none');
    eb2(2) = errorbar(x2(4:6),y2(4:6),yE2(4:6), 'vertical', 'LineStyle', 'none');
    set(eb2, 'color', 'm', 'LineWidth', 1)
    xf = [min(x), max(x)];
    plot(xf,polyval(polyfit(x,y,1), xf),'--k');
    xf2 = [min(x2), max(x2)];
    plot(xf2,polyval(polyfit(x2,y2,1), xf2),'--k');
    ylim([0 3])
    legend('flicker 5hz','flicker 2.5hz','moving 5hz','moving 2.5hz')
    ylabel('motion rating')
    xlabel('SSVEP amplitude')
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'correlSSVEPratingsFqLoc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'correlSSVEPratingsFqLoc' num2str(pickElec(elec)) '.pdf'],'pdf')
    
    corrcoef(x,y)
    [R, P, RL, RU] = corrcoef(x,y);
    % p-values for testing the hypothesis that there is no relationship (if
    % p<0.05 then there is a sig correlation)
    % RL and RU: lower and upper bounds for a 95% confidence interval 
    Rsq = R(1,2).^2;
    disp(['Rsquare static = ' num2str(Rsq)])
    [R, P, RL, RU] = corrcoef(x2,y2);
    Rsq = R(1,2).^2; 
    disp(['Rsquare moving = ' num2str(Rsq)])
    
    
    % figure;hold on;
    % y = [static(2,:) static(3,:)];
    % x = avAmp(6:15);
    % mdl = fitlm(x,y);
    % plot(mdl)
    % y2 = moving([3 5 7 11 15 9 8]);
    % x2 = avAmp(16:22);
    % mdl2 = fitlm(x2,y2);
    % plot(mdl2)
    % ylim([0 3])
    % legend('flickering','','moving')
    % ylabel('motion rating')
    % xlabel('SSVEP amplitude')
    
    
    %%%%%
    % fitlm(x,y)
    % fitlm(x2,y2)
    
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                    Estimate      SE        tStat       pValue
    %                    ________    _______    _______    __________
    %     (Intercept)      2.0997    0.18333     11.453    3.0559e-06
    %     x1             -0.43771    0.11941    -3.6655     0.0063505
    %
    % Number of observations: 10, Error degrees of freedom: 8
    % Root Mean Squared Error: 0.306
    % R-squared: 0.627,  Adjusted R-Squared 0.58
    % F-statistic vs. constant model: 13.4, p-value = 0.00635
    % ans =
    % Linear regression model:
    %     y ~ 1 + x1
    %
    % Estimated Coefficients:
    %                    Estimate       SE        tStat       pValue
    %                    ________    ________    _______    __________
    %     (Intercept)      3.2043      0.1346     23.806    2.4359e-06
    %     x1             -0.55049    0.092074    -5.9787     0.0018757
    %
    % Number of observations: 7, Error degrees of freedom: 5
    % Root Mean Squared Error: 0.179
    % R-squared: 0.877,  Adjusted R-Squared 0.853
    % F-statistic vs. constant model: 35.7, p-value = 0.00188
    
%     %%%%%%% CORREL ALL FREQ
%     y = [static(1,:) static(2,:) static(3,:)];
%     x = avAmp(1:15);
%     y2 = moving([7 5 8 11 3 9 15]); 
%     x2 = avAmp([18 17 22 19 16 21 20]);
    
end

