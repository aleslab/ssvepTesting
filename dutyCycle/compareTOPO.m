load dataDCavRef.mat
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

% plot topo templates
load biosemi128.mat % load biosemi128noInterp.mat
loc = [1:9;10:18]; loc = loc(:);
mm = round(max(max(abs(biosemi128.weights))),-1);
figure('position', [200, 1000, 2000, 500])
for roi=1:18
    subplot(2,9,loc(roi))
    title(biosemi128.listROIs(roi))
    plotTopo(biosemi128.weights(:,roi),cfg.layout)
    caxis([-mm mm])
    colorcet('D1')
end
saveas(gcf,'figures/templatesBiosemi128','png')

        
% Need to normalise amplitudes
% load the amplitudes for a square function
load('squareFFTnew.mat') % these are amplitudes (not power)

% normalisation cannot be done on each fq amp: squareFFT is sometimes 0 in
% harmonics (whole point of doing this). Similarly cannot do odd/even
% separately
% Only the sum can be normalised or F1 alone. 
% + normalisation on power is different than on the sum when there is more than 1 fq!! 
% as shown in the next 2 lines:
cond = 2;
% multiple harmonics
aa = sum(avData(cond).amp(:,avData(cond).harmIdx),2)  / sum(sqAnalytic(avData(cond).harmIdx,cond));
bb= sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2)  / sum(sqAnalytic(avData(cond).harmIdx,cond).^2));
% only f1
cc = sum(avData(cond).amp(:,avData(cond).harmIdx(1)),2)  / sum(sqAnalytic(avData(cond).harmIdx(1),cond));
dd= sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx(1)).^2,2)  / sum(sqAnalytic(avData(cond).harmIdx(1),cond).^2));

for cond=1:22
    topoAmp(:,cond,1) = sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2));
    topoAmp(:,cond,2) = sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2)  / sum(sqAnalytic(avData(cond).harmIdx,cond).^2));
    topoAmp(:,cond,3) = avData(cond).amp(:,(avData(cond).harmIdx(1))) ;
    topoAmp(:,cond,4) = avData(cond).amp(:,(avData(cond).harmIdx(1))) ./ sqAnalytic(avData(cond).harmIdx(1),cond);
    topoAmp(:,cond,5) = sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx(2:end)).^2,2));
    topoAmp(:,cond,6) = sqrt(sum(avData(cond).amp(:,avData(cond).harmIdx(2:end)).^2,2)  / sum(sqAnalytic(avData(cond).harmIdx(2:end),cond).^2));    
end
topoCond = {'amp','ampNorm','f1','f1norm','ampNoF1','ampNormNoF1'};
for ntopo=1:6
    figure; hold on;
    colormap('hot')
    for cond=1:15
        subplot(3,5,cond); hold on;
        plotTopo(topoAmp(:,cond,ntopo),cfg.layout);
        colorbar
        if cond < 6
            caxis([0 max(max(topoAmp(:,[1:5 18],ntopo)))]);
        elseif cond > 5 && cond < 11
            caxis([0 max(max(topoAmp(:,[6:10 17 22 19],ntopo)))]);
        elseif cond > 10
            caxis([0 max(max(topoAmp(:,[11:15 16 21 22],ntopo)))]);
        end
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1500 1000])
    saveas(gcf,['figures' filesep 'topoNew' topoCond{ntopo} '.png'])
end

% motion
position = [11 7 3 9 15 13 8];
for ntopo=1:6
    figure; hold on;
    colormap('hot')
    for cond=16:length(avData)
        subplot(3,5,position(cond-15)); hold on;
        plotTopo(topoAmp(:,cond,ntopo),cfg.layout);
        colorbar
        if position(cond-15) < 6
            caxis([0 max(max(topoAmp(:,[1:5 18],ntopo)))]);
        elseif position(cond-15) > 5 && position(cond-15) < 11
            caxis([0 max(max(topoAmp(:,[6:10 17 22 19],ntopo)))]);
        elseif position(cond-15) > 10
            caxis([0 max(max(topoAmp(:,[11:15 16 21 22],ntopo)))]);
        end
    end
    set(gcf, 'Position', [0 0 1500 1000])
    saveas(gcf,['figures' filesep 'topoNewMotion' topoCond{ntopo} '.png'])
end



% compare for fitting templates
% pool cond sin and cos reg parameter SAME for all conditions
poolDataF1 = []; poolDataNoF1=[];poolDataAll=[];
for cond = 1:length(avData)
    poolDataF1 = [poolDataF1 avData(cond).cos(:,avData(cond).harmIdx(1)) avData(cond).sin(:,avData(cond).harmIdx(1))];
    poolDataNoF1 = [poolDataNoF1 avData(cond).cos(:,avData(cond).harmIdx(2:end)) avData(cond).sin(:,avData(cond).harmIdx(2:end))];
    poolDataAll = [poolDataAll avData(cond).cos(:,avData(cond).harmIdx) avData(cond).sin(:,avData(cond).harmIdx)];
end
[betaF1, betaUnscaledF1,lambdaF1] = fitEEGTemplates(poolDataF1,biosemi128,1);
[betaNoF1, betaUnscaledNoF1,lambdaNoF1] = fitEEGTemplates(poolDataNoF1,biosemi128,1);
[betaAll, betaUnscaledAll,lambdaAll] = fitEEGTemplates(poolDataAll,biosemi128,1);

% cos and sin are stacked across conditions
% need to get their indexes for computing power per condition
nbHarm = arrayfun(@(x) numel([avData(x).harmIdx]),1:22);
% get the range: cumulative sum of elements
rangeHarm = [0 cumsum(nbHarm*2)]; % *2 to account for sin and cos

nbHarmF1 = arrayfun(@(x) numel([avData(x).harmIdx(1)]),1:22);
rangeHarmF1 = [0 cumsum(nbHarmF1*2)];
nbHarmNoF1 = arrayfun(@(x) numel([avData(x).harmIdx(2:end)]),1:22);
rangeHarmNoF1 = [0 cumsum(nbHarmNoF1*2)];

% compute sum of power
for cond=1:length(avData)
    sumSqNoF1(cond) = sum(sqAnalytic(avData(cond).harmIdx(2:end),cond).^2);
    sumSq(cond) = sum(sqAnalytic(avData(cond).harmIdx,cond).^2);
end

% compute amplitude (power=cos^2+sin^2)
for cond = 1:22
    power(:,cond,1) = sum(betaUnscaledF1(:,rangeHarmF1(cond)+1:2:rangeHarmF1(cond+1)).^2 ...
        + betaUnscaledF1(:,rangeHarmF1(cond)+2:2:rangeHarmF1(cond+1)).^2,2);% odd betas = cos, even betas = sin
    power(:,cond,2) = sum(betaUnscaledNoF1(:,rangeHarmNoF1(cond)+1:2:rangeHarmNoF1(cond+1)).^2 ...
        + betaUnscaledNoF1(:,rangeHarmNoF1(cond)+2:2:rangeHarmNoF1(cond+1)).^2,2);% odd betas = cos, even betas = sin
    power(:,cond,3) = sum(betaUnscaledAll(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + betaUnscaledAll(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2);% odd betas = cos, even betas = sin

    %  amplitudes
    amp(:,cond,1) = sqrt( power(:,cond,1) );
    amp(:,cond,2) = sqrt( power(:,cond,2) );
    amp(:,cond,3) = sqrt( power(:,cond,3) );
    %  norm amplitudes
    amp(:,cond,4) = sqrt( power(:,cond,1) / sqAnalytic(avData(cond).harmIdx(1),cond).^2);
    amp(:,cond,5) = sqrt( power(:,cond,2) / sumSqNoF1(cond));
    amp(:,cond,6) = sqrt( power(:,cond,3) / sumSq(cond));
end


topoCond = {'F1','NoF1','All','normF1','normNoF1','normAll'};

col={'b','r','g'};
dcVal = [12.5 25 50 75 87.5];
for hcond=1:6
figure;set(gcf, 'Position', [0 0 1600 700])
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,amp(roi,1:5,hcond),['.-'  col{1}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,amp(roi,6:10,hcond),['.-'  col{2}],'MarkerSize',40,'LineWidth',2)
    plot(dcVal,amp(roi,11:15,hcond),['.-'  col{3}],'MarkerSize',40,'LineWidth',2)
    plot(50,amp(roi,18,hcond),['^:' col{1}],'MarkerSize',15,'Linewidth',2)
    plot([25 50 75],amp(roi,[17 22 19],hcond),['^:' col{2}],'MarkerSize',15,'Linewidth',2)
    plot([12.5 50 87.5],amp(roi,[16 21 20],hcond),['^:' col{3}],'MarkerSize',15,'Linewidth',2)    
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
%     ylim([min(amp(:)) max(amp(:))]);
    ylim([min(min(amp(:,:,hcond))) max(max(amp(:,:,hcond)))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
legend('10','5','2.5','10mov','5mov','2.5mov','Location','Best')
saveas(gcf,['figures' filesep 'dcAmp' topoCond{hcond} '.png'])
end

% compare F1, Except F1, all on the same plot
perFq = [1:5; 6:10; 11:15]; nn=[10 5 2];
for ff=1:3
figure;set(gcf, 'Position', [0 0 1600 700])
for hcond=4:6
for roi = 1:18
    subplot(3,6,roi); hold on;
    plot(dcVal,amp(roi,perFq(ff,:),hcond),'.-','MarkerSize',40,'LineWidth',2)
    xticks([0:12.5:100]);
    title(biosemi128.listROIs(roi))
    ampLim = amp(:,perFq(ff,:),4:6);
    ylim([min(ampLim(:)) max(ampLim(:))]);
    ylabel('amplitude (a.u.)')
end
xlabel('Duty Cycle')
end
legend('F1','noF1','all','Location','Best')
saveas(gcf,['figures' filesep 'dcAmpPerFq' num2str(nn(ff)) '.png'])
end

% check that I retrieve the right topographies
for cond=1:22
    retrieve = biosemi128.weights * betaUnscaledF1;
    retrieveAmp(:,cond,1) = sqrt(sum(retrieve(:,rangeHarmF1(cond)+1:2:rangeHarmF1(cond+1)).^2 ...
        + retrieve(:,rangeHarmF1(cond)+2:2:rangeHarmF1(cond+1)).^2,2));% odd betas = cos, even betas = sin
    retrieveAmp(:,cond,4) = sqrt(sum(retrieve(:,rangeHarmF1(cond)+1:2:rangeHarmF1(cond+1)).^2 ...
        + retrieve(:,rangeHarmF1(cond)+2:2:rangeHarmF1(cond+1)).^2,2) / sum(sqAnalytic(avData(cond).harmIdx(1),cond).^2) );% odd betas = cos, even betas = sin

    retrieve = biosemi128.weights * betaUnscaledNoF1;
    retrieveAmp(:,cond,2) = sqrt(sum(retrieve(:,rangeHarmNoF1(cond)+1:2:rangeHarmNoF1(cond+1)).^2 ...
        + retrieve(:,rangeHarmNoF1(cond)+2:2:rangeHarmNoF1(cond+1)).^2,2));% odd betas = cos, even betas = sin
    retrieveAmp(:,cond,5) = sqrt(sum(retrieve(:,rangeHarmNoF1(cond)+1:2:rangeHarmNoF1(cond+1)).^2 ...
        + retrieve(:,rangeHarmNoF1(cond)+2:2:rangeHarmNoF1(cond+1)).^2,2) / sumSqNoF1(cond) );% odd betas = cos, even betas = sin

    retrieve = biosemi128.weights * betaUnscaledAll;
    retrieveAmp(:,cond,3) = sqrt(sum(retrieve(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + retrieve(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2));% odd betas = cos, even betas = sin
    retrieveAmp(:,cond,6) = sqrt(sum(retrieve(:,rangeHarm(cond)+1:2:rangeHarm(cond+1)).^2 ...
        + retrieve(:,rangeHarm(cond)+2:2:rangeHarm(cond+1)).^2,2) / sumSq(cond) );% odd betas = cos, even betas = sin
end

for ntopo=1:6
    figure; hold on;
    colormap('hot')
    for cond=1:15
        subplot(3,5,cond); hold on;
        plotTopo(retrieveAmp(:,cond,ntopo),cfg.layout);
        colorbar
        if cond < 6
            caxis([0 max(max(retrieveAmp(:,[1:5 18],ntopo)))]);
        elseif cond > 5 && cond < 11
            caxis([0 max(max(retrieveAmp(:,[6:10 17 22 19],ntopo)))]);
        elseif cond > 10
            caxis([0 max(max(retrieveAmp(:,[11:15 16 21 22],ntopo)))]);
        end
    end
    set(gcf, 'Position', [0 0 1500 1000])
    saveas(gcf,['figures' filesep 'dcAmpRetrieved' topoCond{ntopo} '.png'])
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Correl RATINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ratings10.mat')
for fq=1:3
    for dc=1:5
        statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(length(tabStatic));
        movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(length(tabMot));
    end
end

dataToCorrel = topoAmp(:,:,4); % F1 norm
dataToCorrel = topoAmp(:,:,6); % noF1 norm

figure;hold on;
scatter(dataToCorrel(chan,1:5),static(1,:),80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,6:10),static(2,:),80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,11:15),static(3,:),80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,[17 22 19]),moving([5 8 11]),80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,[16 21 20]),moving([3 9 15]),80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,22),moving(7),80,'filled','MarkerEdgeColor','none');
legend('10Hz','5Hz','2Hz','5mov','2mov','10mon')
xlabel('norm SSVEP amplitude')
ylim([0 3]);ylabel('motion rating')
saveas(gcf,['figures' filesep 'correlF1norm.png'])

figure;hold on;
scatter(dataToCorrel(chan,6:15), [static(2,:) static(3,:)],80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,[16:17 19:22]), moving([3 5 11 15 9 8]),80,'filled','MarkerEdgeColor','none');
lsline
ylim([0 3])
legend('flickering','moving')
ylabel('motion rating')
xlabel('norm SSVEP amplitude')
R = corrcoef(dataToCorrel(chan,6:15),[static(2,:) static(3,:)]);
RsqS = R(1,2).^2;
R = corrcoef(dataToCorrel(chan,[16:17 19:22]),moving([3 5 11 15 9 8]));
RsqM = R(1,2).^2;
title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
axis square
saveas(gcf,['figures' filesep 'correlRatingsOzF1.png'])

figure;hold on;
scatter(dataToCorrel(chan,1:15), [static(1,:) static(2,:) static(3,:)],80,'filled','MarkerEdgeColor','none');
scatter(dataToCorrel(chan,[16:17 19:22]), moving([3 5 11 15 9 8]),80,'filled','MarkerEdgeColor','none');
lsline
ylim([0 3])
legend('flickering','moving')
ylabel('motion rating')
xlabel('norm SSVEP amplitude')
R = corrcoef(dataToCorrel(chan,6:15),[static(2,:) static(3,:)]);
RsqS = R(1,2).^2;
R = corrcoef(dataToCorrel(chan,[16:17 19:22]),moving([3 5 11 15 9 8]));
RsqM = R(1,2).^2;
title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
axis square

% x=dataToCorrel(chan,6:15); y=[static(2,:) static(3,:)];
% x2=dataToCorrel(chan,[16:17 19:22]); y2=moving([3 5 11 15 9 8]);
% figure;
% subplot(1,3,1);hold on;
% scatter(x, y,80,'filled','MarkerEdgeColor','none');
% scatter(x2, y2,80,'filled','MarkerEdgeColor','none');
% R = corrcoef(x,y); Rm = corrcoef(x2,y2);
% RsqS = R(1,2).^2; RsqM = Rm(1,2).^2;
% title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
% lsline;ylabel('(motion rating)'); xlabel('(norm SSVEP amplitude)')
% legend('flickering','moving')
% subplot(1,3,2);hold on;
% scatter(log(x),y,80,'filled','MarkerEdgeColor','none');
% scatter(log(x2),y2, 80,'filled','MarkerEdgeColor','none');
% R = corrcoef(log(x),y); Rm = corrcoef(log(x2),y2);
% RsqS = R(1,2).^2; RsqM = Rm(1,2).^2;
% title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
% lsline
% legend('flickering','moving')
% ylabel('motion rating')
% xlabel('log(norm SSVEP amplitude)')
% subplot(1,3,3);hold on;
% scatter(log(x), log(y),80,'filled','MarkerEdgeColor','none');
% scatter(log(x2), log(y2),80,'filled','MarkerEdgeColor','none');
% R = corrcoef(log(x),log(y)); Rm = corrcoef(log(x2),log(y2));
% RsqS = R(1,2).^2; RsqM = Rm(1,2).^2;
% title(['Oz S=' num2str(RsqS,2) ' M=' num2str(RsqM,2) ]) %only 2 digits
% lsline
% legend('flickering','moving')
% ylabel('log(motion rating)')
% xlabel('log(norm SSVEP amplitude)')
% set(gcf, 'Position', [200 200 800 300])
% saveas(gcf,['figures' filesep 'correlRatingsOzF1.png'])


for fqCond = 1:6
    figure;hold on;
    for roi = 1:18
    subplot(3,6,roi); hold on;
    scatter(amp(roi,6:15,fqCond), [static(2,:) static(3,:)],'filled','MarkerEdgeColor','none');
    scatter(amp(roi,[16:17 19:22],fqCond), moving([3 5 11 15 9 8]),'filled','MarkerEdgeColor','none');
    lsline
    ylim([0 3])
%     legend('flickering','moving')
    ylabel('motion rating')
    R = corrcoef(amp(roi,6:15,fqCond),[static(2,:) static(3,:)]);
    RsqS = R(1,2).^2;
    R = corrcoef(amp(roi,[16:17 19:22],fqCond),moving([3 5 11 15 9 8]));
    RsqM = R(1,2).^2;
    title([biosemi128.listROIs{roi} ' S' num2str(RsqS,1) ' M' num2str(RsqM,1) ]) %only 1 digit
%     axis square
    end
    set(gcf, 'Position', [0 0 1200 800])
    saveas(gcf,['figures' filesep 'correlRatings' topoCond{fqCond} '.png'])
end

for fqCond = 1:6
    figure;hold on;
    for roi = 1:18
    subplot(3,6,roi); hold on;
    scatter(log(amp(roi,6:15,fqCond)), [static(2,:) static(3,:)],'filled','MarkerEdgeColor','none');
    scatter(log(amp(roi,[16:17 19:22],fqCond)), moving([3 5 11 15 9 8]),'filled','MarkerEdgeColor','none');
    lsline
    ylim([0 3])
%     legend('flickering','moving')
    ylabel('motion rating')
    R = corrcoef(log(amp(roi,6:15,fqCond)),[static(2,:) static(3,:)]);
    RsqS = R(1,2).^2;
    R = corrcoef(log(amp(roi,[16:17 19:22],fqCond)),moving([3 5 11 15 9 8]));
    RsqM = R(1,2).^2;
    title([biosemi128.listROIs{roi} ' S' num2str(RsqS,1) ' M' num2str(RsqM,1) ]) %only 1 digit
%     axis square
    end
    set(gcf, 'Position', [0 0 1200 800])
    saveas(gcf,['figures' filesep 'correlRatingsLog' topoCond{fqCond} '.png'])
end