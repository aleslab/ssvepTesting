
clearvars;
close all;
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions

for ee=1:2
    load(['allRMSe' num2str(ee)]);
    load(['regCoefE' num2str(ee)]);

    noiseCorr = zeros(length(rmsNoise),11);
    noiseCorr2 = zeros(length(rmsNoise),11); % this is just to have a look with the sqrt results

    for ss=1:length(rmsNoise)
        for motrange=1:2
            noiseCorr(ss,1+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2); % AM + lin
            noiseCorr2(ss,1+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3); % AM + lin
            noiseCorr(ss,2+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2); % AM + spat
            noiseCorr2(ss,2+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,3+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % AM + temp
            noiseCorr2(ss,3+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,4+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % AM + spat + temp
            noiseCorr2(ss,4+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(4);
            % noise from regression uses the coef
            noiseCorr(ss,5+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + coefS(ss,motrange)*rmsNoise(ss,3+5*(motrange-1))^2); % AM + spat
            noiseCorr(ss,6+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + coefT(ss,motrange)*rmsNoise(ss,4+5*(motrange-1))^2); % AM + temp
            noiseCorr(ss,7+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + coefF(ss,motrange,1)*rmsNoise(ss,3+5*(motrange-1))^2 + coefF(ss,motrange,2)*rmsNoise(ss,4+5*(motrange-1))^2);             
        end
    end
    
    % apply normalisation by noise for each participant
    normRMSlr = [rmseNoCoef(:,1:4) ./ noiseCorr(:,1:4) ];
    normRMSsr = [rmseNoCoef(:,5:8) ./ noiseCorr(:,8:11) ];
    
    % take the inverse
    normRMSlr = 1 ./ normRMSlr;
    normRMSsr = 1 ./ normRMSsr;
    
    
    figure; hold on;
    subplot(2,1,1);hold on;
    boxplot(normRMSlr)
    ylabel('1/NRMSE')
    line([1 4],[1 1],'Color','r','LineWidth',2)
    title(['E' num2str(ee) ' LR'])
    xticklabels({'RMSElin','RMSEspat','RMSEtemp','RMSEs+t'})
    ylim([0 1.5])
    subplot(2,1,2);hold on;
    boxplot(normRMSsr)
    ylabel('1/NRMSE')
    line([1 4],[1 1],'Color','r','LineWidth',2)
    title(['E' num2str(ee) ' SR'])
    xticklabels({'RMSElin','RMSEspat','RMSEtemp','RMSEs+t'})
    ylim([0 1.5])
    saveas(gcf,['figures' filesep 'RMSE' num2str(ee)],'png')
    

    %%%%%% supplementary with coefficients
    % apply normalisation by noise for each participant
    normRMSlrCoef = [rmseCoef(:,1:3) ./ noiseCorr(:,5:7)];
    normRMSsrCoef = [rmseCoef(:,4:6) ./ noiseCorr(:,12:14)];
    % take the inverse
    normRMSlrCoef = 1 ./ normRMSlrCoef;
    normRMSsrCoef = 1 ./ normRMSsrCoef;
    figure; hold on;
    subplot(2,1,1);hold on;
    boxplot(normRMSlrCoef)
    ylabel('1/NRMSE')
    line([1 3],[1 1],'Color','r','LineWidth',2)
    title(['E' num2str(ee) ' LR'])
    xticklabels({'RMSEspat','RMSEtemp','RMSEs+t'})
    ylim([0 1.5])
    subplot(2,1,2);hold on;
    boxplot(normRMSsrCoef)
    ylabel('1/NRMSE')
    line([1 3],[1 1],'Color','r','LineWidth',2)
    title(['E' num2str(ee) ' SR'])
    xticklabels({'RMSEspat','RMSEtemp','RMSEs+t'})
    ylim([0 1.5])
    saveas(gcf,['figures' filesep 'RMSEcoef' num2str(ee)],'png')

    
    
    
    % tests the hypothesis that data in x has a continuous distribution with zero median 
    comparisons = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    for cond=1:length(comparisons)
        [pSR(cond)] = signtest(normRMSsr(:,comparisons(cond,1)),normRMSsr(:,comparisons(cond,2))); 
        [pLR(cond)] = signtest(normRMSlr(:,comparisons(cond,1)),normRMSlr(:,comparisons(cond,2))); 
    end
%     save(['signTestE' num2str(ee) '.mat'],'pSR','pLR')    
    




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % topographies
    clear noiseTopoCorr rmsTopoNoise;
    load(['chanRMSall' num2str(ee)]); 
    
    % noise correction for the differences (RMSE)
    noiseTopoCorr = zeros(size(rmsTopoNoise,1),6,size(rmsTopoNoise,3));
    for ss=1:size(rmsTopoNoise,1)
        for chan=1:length(rmsTopoNoise)
            cond = 1;
            for spat=[1 5]
            for dd=1:3
                noiseTopoCorr(ss,cond,chan) = sqrt(rmsTopoNoise(ss,spat,chan)^2 + rmsTopoNoise(ss,spat+dd,chan)^2); 
                cond = cond+1;
            end
            end
        end
    end
%     figure; plotTopo(squeeze(mean(noiseTopoCorr(:,4,:))),cfg.layout); colorbar

    % apply normalisation
    % rmsTopo in the order: rms AM, rms Lin, rmse Lin, rmse Spat, rmse Temp *2
    % rmsTopoNoise in the order AM, Lin, spat, temp  *2 
    normTopoRMSlr = [rmsTopo(:,1:2,:)./rmsTopoNoise(:,1:2,:) rmsTopo(:,3:5,:) ./ noiseTopoCorr(:,1:3,:) ];
    normTopoRMSsr = [rmsTopo(:,6:7,:)./rmsTopoNoise(:,5:6,:) rmsTopo(:,8:10,:) ./ noiseTopoCorr(:,4:6,:) ];
    
    
    % plot topo
    cfg.layout = 'biosemi128.lay';
    titre = {'RMS AM','RMS Lin','RMSE Lin', 'RMSE spat','RMSE temp'};
    figure('Renderer', 'painters', 'Position', [10 10 1400 700])
    for cond=1:5
        subplot(2,5,cond)
        plotTopo(squeeze(mean(normTopoRMSlr(:,cond,:))),cfg.layout)
        colorbar
        if ee==1
            caxis([1 2.5])
        elseif ee == 2 && cond<3
            caxis([1 7.5])
        else
            caxis([1 3])
        end
        title(titre{cond})
        subplot(2,5,cond+5)
        plotTopo(squeeze(mean(normTopoRMSsr(:,cond,:))),cfg.layout)
        colorbar
        if ee==1
            caxis([1 2.5])
        elseif ee == 2 && cond<3
            caxis([1 7.5])
        else
            caxis([1 3])
        end
        title(titre{cond})
    end
    colormap(jmaColors('hotcortex'));
    colormap('hot');
    saveas(gcf,['figures' filesep 'topoRMS E' num2str(ee)],'png')
    
    
    %%%% RMS Electrodes picked
    pickElec = [23 126 38];

        normTopoRMSlr = [rmsTopo(:,1:2,:)./rmsTopoNoise(:,1:2,:) rmsTopo(:,3:5,:) ./ noiseTopoCorr(:,1:3,:) ];
    normTopoRMSsr = [rmsTopo(:,6:7,:)./rmsTopoNoise(:,5:6,:) rmsTopo(:,8:10,:) ./ noiseTopoCorr(:,4:6,:) ];
    
    figure; hold on;
    for chan=1:length(pickElec)
        subplot(3,1,chan); hold on;
        tmpRMSlr = [rmsTopo(:,1:2,pickElec(chan))./rmsTopoNoise(:,1:2,pickElec(chan)) rmsTopo(:,3:5,pickElec(chan)) ./ noiseTopoCorr(:,1:3,pickElec(chan))];
        tmpRMSsr = [rmsTopo(:,6:7,pickElec(chan))./rmsTopoNoise(:,6:7,pickElec(chan)) rmsTopo(:,8:10,pickElec(chan)) ./ noiseTopoCorr(:,4:6,pickElec(chan))];
        boxplot([1./tmpRMSlr 1./tmpRMSsr])
        line([1 10],[1 1],'Color','r','LineWidth',2)
        title(['elec' num2str(pickElec(chan))])
        ylabel('1/ (RMS S/N)')
        ylim([-0.2 2.2])
        xticklabels({'longRMSAM','RSMlin','RMSElin','RMSEspa','RMSEtmp','shortAM','RSMlin','RMSElin','RMSEspa','RMSEtmp'})
    end
    saveas(gcf,['figures' filesep 'RMS per elec E' num2str(ee)],'png')
    
    
end
