addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% see pred.labels for details
% first 12 are 1st exp, from 12 is the 2nd exp
clearvars;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};



addpath /Users/Marlene/Documents/git/fieldtrip-aleslab-fork
addpath /Users/Marlene/Documents/git/ssvepTesting/svndlCopy
addpath /Users/Marlene/Documents/git/ssvepTesting/biosemiUpdated
addpath /Users/Marlene/Documents/git/ssvepTesting/commonFunctions
ft_defaults

% see pred.labels for details
% first 12 are 1st exp, from 12 is the 2nd exp
clearvars;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

dataPath = {'/Users/Marlene/Documents/git/dataLR/LRlongDC/', '/Users/Marlene/Documents/git/dataLR/LRshortDC/'};


for ee=1:2 % which experiment
    clear sbj
    load([dataPath{ee} 'sbjprediction.mat'])
    
    % conditions in sbj: AM, linear, spatial, temp, S+T, AM short, linear,
    % spatial, temp, S+T, non-linear ST long, non-linear ST short
    
    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% over all time and electrodes
    clear rmsNoise rmsSignal rmse nrms nrmse nrmseCoef Bnrmse
    % RMS over all time and electrodes
    for ss=1:length(sbj)
        for numCond=1:10 
            rmsNoise(ss,numCond) = rms(sbj(ss,numCond).data.noiseWave(:)  );
            rmsSignal(ss,numCond) = rms(sbj(ss,numCond).data.filteredWave(:)  );
        end
    end
    % RMSE over all time and electrodes
    % lin, spa, temp, s+t
    for ss=1:length(sbj)
        for numCond=2:5 
            rmse(ss,numCond-1) = rms(sbj(ss,1).data.filteredWave(:) -  sbj(ss,numCond).data.filteredWave(:)); % long-range
            rmse(ss,numCond+3) = rms(sbj(ss,6).data.filteredWave(:) -  sbj(ss,numCond+5).data.filteredWave(:)); % short-range
        end
    end
    % NRMS = RMS / RMS noise
    for ss=1:length(sbj)
        for numCond=1:size(rmsSignal,2)
            nrms(ss,numCond) = rmsSignal(ss,numCond) / rmsNoise(ss,numCond) ;
        end
    end    
    % NRMSE = RMSE / RMS noise depending on the condition
    for ss=1:length(sbj)
        for numCond=2:4
            nrmse(ss,numCond-1) = rmse(ss,numCond-1) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,numCond)^2)); % divide by noise AM + noise from the condition
            nrmse(ss,numCond+3) = rmse(ss,numCond+3) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,numCond+5)^2)); % divide by noise AM + noise from the condition
        end
        % for s+t need to add both noise + AM
        nrmse(ss,4) = rmse(ss,4) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,3)^2 + rmsNoise(ss,4)^2)); 
        nrmse(ss,8) = rmse(ss,8) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,8)^2 + rmsNoise(ss,9)^2)); 
    end
    % NRMSE after regression for LR and SR (done only for aS+bT)
    load(['rmseCoefE' num2str(ee) '.mat']) % rmseS rmseT rmseS+T x LR/SR 
    load(['regCoefE' num2str(ee)]); % get the regression coef to get the amount of noise included
    for ss=1:length(sbj)
        nrmseCoef(ss,1) = rmseCoef(ss,3) / (sqrt(rmsNoise(ss,1)^2 + coefF(ss,1,1)*rmsNoise(ss,3)^2 + coefF(ss,1,2)*rmsNoise(ss,4)^2)); % CHECK it is coef*(N^2) not (coef*N)^2
        nrmseCoef(ss,2) = rmseCoef(ss,6) / (sqrt(rmsNoise(ss,6)^2 + coefF(ss,2,1)*rmsNoise(ss,8)^2 + coefF(ss,2,2)*rmsNoise(ss,9)^2 )); 
    end
    
%     %%% plot NMRSE
%     figure('Renderer', 'painters', 'Position', [10 10 400 900])
%     hold on;
%     subplot(2,1,1)
%     boxplot(1./[nrmse(:,1:4) nrmseCoef(:,1)])
%     line([1 5],[1 1],'Color','r','LineWidth',2)
%     title('LR')
%     ylim([0 1.2])
%     xticklabels({'lin','spat','temp','s+t','regress'})
%     ylabel('1/NRMSE')
%     subplot(2,1,2)
%     boxplot(1./[nrmse(:,5:8) nrmseCoef(:,2)])
%     title('SR')
%     line([1 5],[1 1],'Color','r','LineWidth',2)
%     xticklabels({'lin','spat','temp','s+t','regress'})
%     ylabel('1/NRMSE')
%     ylim([0 1.2])
%     saveas(gcf,['figures' filesep 'E' num2str(ee) 'NMRSEboxplot.png'])
    
    %%% use gramm to plot nicer boxplot
    clear g;
    X=[1 2 3 4 5];Y=[1 1 1 1 1];
    yval = [nrmse(:,1:4) nrmseCoef(:,1) nrmse(:,5:8) nrmseCoef(:,2)];
    xval = repmat({'1lin';'2spat';'3temp';'4s+t';'5regress'},1,length(nrmse))'; % numbers to avoid alphabetic order
    g(1,1) = gramm('x',xval(:),'y',1./yval(1:end/2),'ymin',repmat(0,length(yval(1:end/2)),1),'ymax',repmat(1.2,length(yval(1:end/2)),1));
    g(2,1) = gramm('x',xval(:),'y',1./yval(end/2+1:end),'ymin',repmat(0,length(yval(1:end/2)),1),'ymax',repmat(1.2,length(yval(1:end/2)),1));
    g(1,1).stat_boxplot();g(2,1).stat_boxplot();
    g(1,1).set_names('x','prediction','y','1/NMRSE');
    g(2,1).set_names('x','prediction','y','1/NMRSE');
    g(2,1).set_title('SR');
    g(1,1).set_title('LR');
    g(1,1).geom_hline('yintercept',1);g(2,1).geom_hline('yintercept',1)
    figure('Renderer', 'painters', 'Position', [10 10 400 600])
    g.draw();
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'NMRSEboxplot.png'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'NMRSEboxplot.pdf'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'NMRSEboxplot.fig'])
    save(['statNRMSE' num2str(ee) '.mat'],'yval') % save NRMSE for stats


%     % instead of boxplot, determine 95% CI
%     % because of the hard limit at 1 the distribution is skewed (F
%     % distribution) so more informative to show 95%CI
%     % for this, get a distrib for 1000 bootstrap of sbj with replacement
%     for boot=1:1000
%         nn=0;
%         for ss=randsample(1:length(rmse),length(rmse),true)
%             nn=nn+1;
%             for numCond=2:4
%                 Bnrmse(nn,numCond-1) = rmse(ss,numCond-1) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,numCond)^2)); % divide by noise AM + noise from the condition
%                 Bnrmse(nn,numCond+3) = rmse(ss,numCond+3) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,numCond+5)^2)); % divide by noise AM + noise from the condition
%             end
%             % for s+t need to add both noise + AM
%             Bnrmse(nn,4) = rmse(ss,4) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,3)^2 + rmsNoise(ss,4)^2));
%             Bnrmse(nn,8) = rmse(ss,8) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,8)^2 + rmsNoise(ss,9)^2));
%             % for coef
%             Bnrmse(nn,9) = rmseCoef(ss,3) / (sqrt(rmsNoise(ss,1)^2 + coefF(ss,1,1)*rmsNoise(ss,3)^2 + coefF(ss,1,2)*rmsNoise(ss,4)^2));
%             Bnrmse(nn,10) = rmseCoef(ss,6) / (sqrt(rmsNoise(ss,6)^2 + coefF(ss,2,1)*rmsNoise(ss,8)^2 + coefF(ss,2,2)*rmsNoise(ss,9)^2 ));
%         end
%         distNRMSE(boot,:) = mean(1./Bnrmse);
%     end
%     for cond=1:size(distNRMSE,2)
%         ciNRMSE(cond,:) = quantile(distNRMSE(:,cond),[0.025 0.975]);
%     end
%     figure('Renderer', 'painters', 'Position', [10 10 900 400])
%     subplot(1,2,1);hold on;
%     for cond=1:4
%     errorbar(cond,mean(1./nrmse(:,cond)),mean(1./nrmse(:,cond))-ciNRMSE(cond,1),ciNRMSE(cond,2)-mean(1./nrmse(:,cond)),'Marker','o','CapSize',0);
%     end
%     ylim([0 1])
%     xlim([0.5 4.5])
%     xticklabels({'lin','spat','temp','s+t'})
%     ylabel('1/NRMSE')
%     title('LR')
%     subplot(1,2,2);hold on;
%     for cond=5:8
%     errorbar(cond,mean(1./nrmse(:,cond)),mean(1./nrmse(:,cond))-ciNRMSE(cond,1),ciNRMSE(cond,2)-mean(1./nrmse(:,cond)),'Marker','o','CapSize',0);
%     end
%     ylim([0 1])
%     xlim([4.5 8.5])
%     ylabel('1/NRMSE')
%     title('SR')
%     xticklabels({'lin','spat','temp','s+t'})
%     saveas(gcf,['E' num2str(ee) 'NMRSE95CI.png'])
  


    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Per electrode
    clear rmsNoiseChan rmsSignalChan rmseChan nrmseChan nrmseChanCoef
    % RMS 
    for ss=1:length(sbj)
        for numCond=1:10
            for chan = 1 : size(sbj(ss,numCond).data.filteredWave,1)
                rmsNoiseChan(ss,numCond,chan) = rms( sbj(ss,numCond).data.noiseWave(chan,:) );
                rmsSignalChan(ss,numCond,chan) = rms( sbj(ss,numCond).data.filteredWave(chan,:) );
            end
        end
    end
    % RMSE 
    % spa, temp, s+t
    for ss=1:length(sbj)
        for numCond=2:5
            for chan = 1 : size(sbj(ss,numCond).data.filteredWave,1)
                rmseChan(ss,numCond-1,chan) = rms(sbj(ss,1).data.filteredWave(chan,:) -  sbj(ss,numCond).data.filteredWave(chan,:)); % long-range
                rmseChan(ss,numCond+3,chan) = rms(sbj(ss,6).data.filteredWave(chan,:) -  sbj(ss,numCond+5).data.filteredWave(chan,:)); % short-range
            end
        end
    end
    
    %%% plot topographies RMS AM
    figure('Renderer', 'painters', 'Position', [10 10 400 900])
    toPlot = [1 6]; % plot only AM condition
    titre = {'LR','SR'};
    for cond=1:length(toPlot)
        subplot(2,1,cond)
        plotTopo(squeeze(mean(rmsSignalChan(:,toPlot(cond),:))),cfg.layout)
%         plotTopo(squeeze(mean(rmsSignalChan(:,toPlot(cond),:)./rmsNoiseChan(:,toPlot(cond),:))),cfg.layout)
        colorbar
        colormap('hot');
        caxis([0 2.3]) % max is 2.26 in shortDC
        title(titre{cond})
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoAM.png'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoAM.pdf'])
    
    
    %%% plot topographies NRMSE
    for cond=2:5
        nrmseChan(:,cond-1,:) = rmseChan(:,cond-1,:) ./ (sqrt(rmsNoiseChan(:,1,:).^2 + rmsNoiseChan(:,cond,:).^2) ); % LR
        nrmseChan(:,cond-1+4,:) = rmseChan(:,cond-1+4,:) ./ (sqrt(rmsNoiseChan(:,6,:).^2 + rmsNoiseChan(:,cond+5,:).^2) ); % SR
    end
    % NRMSE after regression for LR and SR
    load(['rmseFuTopo' num2str(ee) '.mat']) % LR/SR, same coef as for the across electrodes 
    for ss=1:length(sbj)
        nrmseChanCoef(ss,1,:) = rmseFuTopo(ss,1,:) ./ sqrt( rmsNoiseChan(ss,1,:).^2 + coefF(ss,1,1)*rmsNoiseChan(ss,3,:).^2 + coefF(ss,1,2)*rmsNoiseChan(ss,4,:).^2 ) ; % LR
        nrmseChanCoef(ss,2,:) = rmseFuTopo(ss,2,:) ./ sqrt( rmsNoiseChan(ss,6,:).^2 + coefF(ss,1,1)*rmsNoiseChan(ss,7,:).^2 + coefF(ss,1,2)*rmsNoiseChan(ss,8,:).^2 ) ; % LR
    end

    allNRMSEtopo = [nrmseChan(:,1:4,:) nrmseChanCoef(:,1,:) nrmseChan(:,5:8,:) nrmseChanCoef(:,2,:)];
    figure('Renderer', 'painters', 'Position', [10 10 1600 800])
    titre = {'lin','spat','temp','s+t','regress','lin','spat','temp','s+t','regress'};
    for cond=1:size(allNRMSEtopo,2)
        subplot(2,5,cond)
        plotTopo(squeeze(mean(allNRMSEtopo(:,cond,:))),cfg.layout)
        colorbar
        colormap('hot');
        caxis([1 3.4]) % max is 3.39 in shortDC, min cannot be less than 1
        title(titre{cond})
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoNRMSE.png'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoNRMSE.fig'])
    
    allRMSEtopo = [rmseChan(:,1:4,:) rmseFuTopo(:,1,:) rmseChan(:,5:8,:) rmseFuTopo(:,2,:)];
%     figure('Renderer', 'painters', 'Position', [10 10 1000 500])
figure;
    titre = {'lin','spat','temp','s+t','regress','lin','spat','temp','s+t','regress'};
    for cond=1:size(allRMSEtopo,2)
        subplot(2,5,cond)
        plotTopo(squeeze(mean(allRMSEtopo(:,cond,:))),cfg.layout)
        colorbar
        colormap('hot');
        caxis([0 2.3]) % max is 2.05 in shortDC, min cannot be less than 0
        title(titre{cond})
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoRMSE.png'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'topoRMSE.fig'])
end

%%% stats
cond1 = [1 1 1 1 2 2 2 3 3 4 ];
cond2 = [2 3 4 5 3 4 5 4 5 5 ];
for cc=1:length(cond1)
    pSign(cc) = signtest(1./yval(:,cond1(cc)),1./yval(:,cond2(cc)));
    pSign(cc+ length(cond1)) = signtest(1./yval(:,cond1(cc)+ 5),1./yval(:,cond2(cc)+ 5));
end
