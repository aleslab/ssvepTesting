
clearvars;
close all;
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions

for ee=1:2
    load(['allRMSe' num2str(ee)]);
    load(['regCoefE' num2str(ee)]);
    figure; hold on
    subplot(2,2,1);hold on;
    bar([mean(rmsSignal(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))])
    errorbar([mean(rmsSignal(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))], ...
        [std(rmsSignal(:,1:2)) std(rmseNoCoef(:,1:4)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['rmsE' num2str(ee) 'LR'])
    subplot(2,2,2);hold on;
    bar([mean(rmsSignal(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))])
    errorbar([mean(rmsSignal(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))], ...
        [std(rmsSignal(:,6:7)) std(rmseNoCoef(:,5:8)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['rmsE' num2str(ee) 'SR'])
    % x-axis: AM lin rmseNoCoef (lin s t s+t) rmseCoef (s t s+t)
    
    subplot(2,2,3);hold on;
    bar(mean(rmsNoise))
    errorbar(mean(rmsNoise), std(rmsNoise),'k','Linestyle','none')
    title(['raw noiseE' num2str(ee)])
    ylabel('rms noise per cond')
    % x-axis: LR for AM lin S T ST, SR for AM lin S T ST
    
    % correct noise
    % order: linear spatial temporal s+t
    % eg noiseLin = sqrt(noiseAM*noiseAM + noiseLin*noiseLin)
    % should be similar to noiseAM * sqrt(3)
    % the std per participant for the noise is around 0.5 calculated by:
    % sqrt([(noisecorrected / sqrt(2/pi))^2] * (4-pi)/2)
    noiseCorr = zeros(length(rmsNoise),11);
    noiseCorr2 = zeros(length(rmsNoise),11); % this is just to have a look with the sqrt results

    %%% not including AM noise for the rmse
    for ss=1:length(rmsNoise)
        for motrange=1:2
            noiseCorr(ss,1+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2); % AM + lin
            noiseCorr2(ss,1+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3); % AM + lin
            noiseCorr(ss,2+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2); % lin + spat
            noiseCorr2(ss,2+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,3+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % lin + temp
            noiseCorr2(ss,3+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,4+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % lin + spat + temp
            noiseCorr2(ss,4+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(4);
            % noise from regression uses the coef
            noiseCorr(ss,5+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefS(ss,motrange)*rmsNoise(ss,3+5*(motrange-1))^2); % lin + spat
            noiseCorr(ss,6+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefT(ss,motrange)*rmsNoise(ss,4+5*(motrange-1))^2); % lin + temp
            noiseCorr(ss,7+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefF(ss,motrange,1)*rmsNoise(ss,3+5*(motrange-1))^2 + coefF(ss,motrange,2)*rmsNoise(ss,4+5*(motrange-1))^2);             
        end
    end
    subplot(2,2,4);hold on;
    bar(mean(noiseCorr))
    errorbar(mean(noiseCorr), std(noiseCorr),'k','Linestyle','none')
    title(['corrected noise no coef E' num2str(ee)])
    ylabel('rms noise')
    
    % apply normalisation by noise for each participant
    normRMSlr = [rmsSignal(:,1:2)./rmsNoise(:,1:2) rmseNoCoef(:,1:4) ./ noiseCorr(:,1:4) rmseCoef(:,1:3) ./ noiseCorr(:,5:7)];
    normRMSsr = [rmsSignal(:,6:7)./rmsNoise(:,6:7) rmseNoCoef(:,5:8) ./ noiseCorr(:,8:11) rmseCoef(:,4:6) ./ noiseCorr(:,12:14)];
    
    % take the inverse
    normRMSlr = 1 ./ normRMSlr;
    normRMSsr = 1 ./ normRMSsr;
    
%     figure; hold on;
%     subplot(2,1,1);hold on;
%     bar(mean(normRMSlr))
%     errorbar(mean(normRMSlr),std(normRMSlr)/sqrt(length(normRMSsr)),'k','Linestyle','none')
%     ylabel('RMS signal/ rms noise')
%     line([1 9],[1 1],'Color','r','LineWidth',3)
%     title(['E' num2str(ee) ' LR'])
%     subplot(2,1,2);hold on;
%     bar(mean(normRMSsr))
%     errorbar(mean(normRMSsr),std(normRMSsr)/sqrt(length(normRMSsr)),'k','Linestyle','none')
%     ylabel('RMS signal/ rms noise')
%     line([1 9],[1 1],'Color','r','LineWidth',3)
%     title(['E' num2str(ee) ' SR'])
%     saveas(gcf,['figures' filesep 'RMScorrected' num2str(ee)],'png')
    
    figure; hold on;
    subplot(2,1,1);hold on;
    boxplot(normRMSlr)
    ylabel('1/ (RMS S/N)')
    line([1 9],[1 1],'Color','r','LineWidth',2)
    title(['E' num2str(ee) ' LR'])
    xticklabels({'RMSAM','RSMlin','RMSElin','RMSEspat','RMSEtemp','RMSEst','coefspat','coeftemp','coefST'})
    ylim([0 1.5])
    subplot(2,1,2);hold on;
    boxplot(normRMSsr)
    ylabel('1/ (RMS S/N)')
    line([1 9],[1 1],'Color','r','LineWidth',3)
    title(['E' num2str(ee) ' SR'])
    xticklabels({'RMSAM','RSMlin','RMSElin','RMSEspat','RMSEtemp','RMSEst','coefspat','coeftemp','coefST'})
    ylim([0 1.5])
    saveas(gcf,['figures' filesep 'RMScorrected' num2str(ee)],'png')
    
    
    % tests the hypothesis that data in x has a continuous distribution with zero median 
    for cond=1:size(normRMSsr,2)
        [pSR(cond)] = signtest(normRMSsr(:,cond)-1); 
        [pLR(cond)] = signtest(normRMSlr(:,cond)-1); 
    end
    save(['signTestE' num2str(ee) '.mat'],'pSR','pLR')    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % topographies
    clear noiseTopoCorr rmsTopoNoise;
    load(['chanRMS' num2str(ee)]); 
    load(['rmsFuTopo' num2str(ee)]); 
    
    % noise correction for the differences (RMSE)
    noiseTopoCorr = zeros(size(rmsTopoNoise,1),6,size(rmsTopoNoise,3));
    for ss=1:size(rmsTopoNoise,1)
        for chan=1:length(rmsTopoNoise)
            noiseTopoCorr(ss,1,chan) = sqrt(rmsTopoNoise(ss,1,chan)^2 + rmsTopoNoise(ss,2,chan)^2); % AM + lin
            noiseTopoCorr(ss,2,chan) = sqrt(rmsTopoNoise(ss,2,chan)^2 + rmsTopoNoise(ss,4,chan)^2); % lin + temp
            noiseTopoCorr(ss,3,chan) = sqrt(rmsTopoNoise(ss,5,chan)^2 + rmsTopoNoise(ss,6,chan)^2); % AM + lin
            noiseTopoCorr(ss,4,chan) = sqrt(rmsTopoNoise(ss,6,chan)^2 + rmsTopoNoise(ss,8,chan)^2); % lin + temp
            
            % regression coef 
            noiseTopoCorr(ss,5,chan) = sqrt(rmsTopoNoise(ss,2,chan)^2 + coefF(ss,1,1)*rmsTopoNoise(ss,3,chan)^2 + coefF(ss,1,2)*rmsTopoNoise(ss,4,chan)^2); % lin + spat + temp             
            noiseTopoCorr(ss,6,chan) = sqrt(rmsTopoNoise(ss,6,chan)^2 + coefF(ss,2,1)*rmsTopoNoise(ss,7,chan)^2 + coefF(ss,2,2)*rmsTopoNoise(ss,8,chan)^2); % lin + spat + temp             
        end
    end
%     figure; plotTopo(squeeze(mean(noiseTopoCorr(:,4,:))),cfg.layout); colorbar

    % apply normalisation
    % rmsTopo in the order: rms AM, rms Lin, rmse Lin, rmse Temp
    normTopoRMSlr = [rmsTopo(:,1:2,:)./rmsTopoNoise(:,1:2,:) rmsTopo(:,3:4,:) ./ noiseTopoCorr(:,1:2,:) rmsFuTopo(:,1,:)./noiseTopoCorr(:,5,:)];
    normTopoRMSsr = [rmsTopo(:,5:6,:)./rmsTopoNoise(:,5:6,:) rmsTopo(:,7:8,:) ./ noiseTopoCorr(:,3:4,:) rmsFuTopo(:,2,:)./noiseTopoCorr(:,6,:)];
    
    
    % plot topo
    cfg.layout = 'biosemi128.lay';
    titre = {'RMS AM','RMS Lin','RMSE Lin', 'RMSE temp','RMSE coef'};
    figure('Renderer', 'painters', 'Position', [10 10 1400 700])
    for cond=1:5
        subplot(2,5,cond)
        plotTopo(squeeze(mean(normTopoRMSlr(:,cond,:))),cfg.layout)
        colorbar
        caxis([1 7.5])
        if ee==2
            caxis([1 7.5])
        else
            caxis([1 2.5])
        end
        title(titre{cond})
        subplot(2,5,cond+5)
        plotTopo(squeeze(mean(normTopoRMSsr(:,cond,:))),cfg.layout)
        colorbar
        if ee==2
            caxis([1 7.5])
        else
            caxis([1 2.5])
        end
        title(titre{cond})
    end
    colormap(jmaColors('hotcortex'));
%     colormap('hot');
    saveas(gcf,['figures' filesep 'topoRMS E' num2str(ee)],'png')
    
    
    %%%% RMS Electrodes picked
    pickElec = [23 126 38];

    figure; hold on;
    for chan=1:length(pickElec)
        subplot(3,1,chan); hold on;
        tmpRMSlr = [rmsTopo(:,1:2,pickElec(chan))./rmsTopoNoise(:,1:2,pickElec(chan)) rmsTopo(:,3:4,pickElec(chan)) ./ noiseTopoCorr(:,1:2,pickElec(chan))];
        tmpRMSsr = [rmsTopo(:,5:6,pickElec(chan))./rmsTopoNoise(:,5:6,pickElec(chan)) rmsTopo(:,7:8,pickElec(chan)) ./ noiseTopoCorr(:,3:4,pickElec(chan))];
        boxplot([1./tmpRMSlr 1./tmpRMSsr])
        line([1 8],[1 1],'Color','r','LineWidth',2)
        title(['elec' num2str(pickElec(chan))])
        ylabel('1/ (RMS S/N)')
        ylim([-0.2 2.2])
        xticklabels({'longRMSAM','RSMlin','RMSElin','RMSEtmp','shortAM','RSMlin','RMSElin','RMSEtmp'})
    end
    saveas(gcf,['figures' filesep 'RMS per elec E' num2str(ee)],'png')
    
%     figure; hold on;
%     for chan=1:length(pickElec)
%         subplot(3,1,chan); hold on;
%         tmpRMSlr = [rmsTopo(:,1:2,pickElec(chan))./rmsTopoNoise(:,1:2,pickElec(chan)) rmsTopo(:,3:4,pickElec(chan)) ./ noiseTopoCorr(:,1:2,pickElec(chan))];
%         tmpRMSsr = [rmsTopo(:,5:6,pickElec(chan))./rmsTopoNoise(:,5:6,pickElec(chan)) rmsTopo(:,7:8,pickElec(chan)) ./ noiseTopoCorr(:,3:4,pickElec(chan))];
%         bar([mean(tmpRMSlr) mean(tmpRMSsr)])
%         errorbar([mean(tmpRMSlr) mean(tmpRMSsr)],[std(tmpRMSlr)/sqrt(length(tmpRMSlr)) std(tmpRMSsr)/sqrt(length(tmpRMSlr))],'k','Linestyle','none')
%         line([1 8],[1 1],'Color','r','LineWidth',3)
%         title(['elec' num2str(pickElec(chan))])
%     end
%     saveas(gcf,['figures' filesep 'RMS per elec E' num2str(ee)],'png')
    
    
end
