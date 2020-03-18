
clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% load individual data
dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7:11 13:15 17:18 20]; 
% reject S7 and S13 S17, S20 same as S03
% S1=pilote so reject numS-1 (6 12 16 19)

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    dataSbj(:,ss) = Axx;
end

[avData, proj_Amp] = averageAxxWithStd(dataSbj);

col={'b','r','g'};

% pickElec = 23; % 23 Oz, 9, B6=38, 16=best SNR

% from topo max: A15 & A27
pickElec = [15 23 27];

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Topography
% flicker
figure; hold on;
colormap('hot')
for cond=1:15
%     max(avData(cond).amp(:,avData(cond).i1f1))
    subplot(3,5,cond); hold on;
    max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    plotTopo(avData(cond).amp(:,avData(cond).i1f1),cfg.layout);
    colorbar
    if cond < 6
        caxis([0 1.5]);
    else
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoFlickF1.jpg'])
% motion
position = [11 7 3 9 15 13 8];
figure; hold on;
colormap('hot')
for cond=16:length(avData)
%     max(avData(cond).amp(:,avData(cond).i1f1*2-1))
    subplot(3,5,position(cond-15)); hold on;
    plotTopo(avData(cond).amp(:,avData(cond).i1f1*2-1),cfg.layout);
    colorbar
    if cond-15 == 3
        caxis([0 1.5]);
    else
        caxis([0 3]);
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 1500 1000])
saveas(gcf,['figures' filesep 'topoMotionF2.jpg'])

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SSVEP

% get F1 amplitude
for elec = 1:length(pickElec)
    for cond=1:15
        avAmp(cond)=avData(cond).amp(pickElec(elec),(avData(cond).i1f1));
        %    avProj(cond)=mean(proj_Amp(pickElec,avData(cond).i1f1,cond,:));
        semProj(cond) = std(proj_Amp(pickElec(elec),avData(cond).i1f1,cond,:)) / sqrt(size(dataSbj,2));
    end
    for cond=16:22 % for motion take the 2f1
        avAmp(cond)=avData(cond).amp(pickElec(elec),(avData(cond).i1f1*2-1));
        semProj(cond) = std(proj_Amp(pickElec(elec),avData(cond).i1f1*2-1,cond,:)) / sqrt(size(dataSbj,2));
    end
    % % compute noise level per freq
    % for cond=1:15 % exclude motion conditions
    %    baseline(cond) = (avData(cond).amp(pickElec,avData(cond).i1f1-1)+avData(cond).amp(pickElec,avData(cond).i1f1+1)) / 2;
    % end
    % compute noise level all cond pooled
    for cond = 1:length(avData)
        baseline(cond) = (avData(cond).amp(pickElec(elec),avData(cond).i1f1-1)+avData(cond).amp(pickElec(elec),avData(cond).i1f1+1)) / 2;
    end
    temp=reshape(baseline(1:20),4,5);
    avBaseline = mean(temp);
    avBaseline(3) = mean([temp(1:4,3)' baseline(21) baseline(22)]); % add the 50% DC for other freq
    
    dcVal = [12.5 25 50 75 87.5];
    
    figure;hold on;
    for freq=1:3
        errorbar(dcVal,avAmp((freq-1)*5+1:freq*5),semProj((freq-1)*5+1:freq*5),['.-' col{freq}],'MarkerSize',40,'Linewidth',2);
    end
    plot(dcVal,avBaseline,'--k','Linewidth',2)
    errorbar(50,avAmp(18),semProj(18),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    errorbar([25 50 75],avAmp([17 22 19]),semProj([17 22 19]),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    errorbar([12.5 50 87.5],avAmp([16 21 20]),semProj([17 22 19]),['^:' col{3}],'MarkerSize',15,'Linewidth',2)
    xlim([0 100]);
    % xticks([1:1:5]);
    % xticklabels({'12.5','25','50','75','87.5'})
    ylim([0 inf])
    legend('10','5','2.5','noise level','Location','Best')
    xlabel('Duty Cycle')
    ylabel('SSVEP amplitude')
    set(gca,'FontSize',15)
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'ampDCloc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'ampDCloc' num2str(pickElec(elec)) '.eps'],'epsc')
    
    %%%%%% LINEAR REGRESSION
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
    load('fullRatings.mat')
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
        % xticks([1:1:5]);
        ylim([0 inf])
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
    x = avAmp(6:15);
    y2 = moving([3 5 11 15 9 8]); % do not include 10Hz
    x2 = avAmp([16:17 19:22]);% do not include 10Hz
    
    
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
    scatter(x,y,80,'filled','MarkerEdgeColor','none')
    scatter(x2,y2,80,'filled','MarkerEdgeColor','none')
    xf = [min(x), max(x)];
    plot(xf,polyval(polyfit(x,y,1), xf),'--k');
    xf2 = [min(x2), max(x2)];
    plot(xf2,polyval(polyfit(x2,y2,1), xf2),'--k');
    ylim([0 3])
    legend('flickering','moving')
    ylabel('motion rating')
    xlabel('SSVEP amplitude')
    title(num2str(pickElec(elec)))
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.png'])
    saveas(gcf,['figures' filesep 'correlSSVEPratingsLoc' num2str(pickElec(elec)) '.pdf'],'pdf')
    
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
    
    
end

