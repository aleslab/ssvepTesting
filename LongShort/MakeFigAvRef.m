% Make figures with average Reference
% predictions (with AM)
% (for RMS see plotRMS)

% addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
% ft_defaults

addpath C:\Users\Marlene\Documents\git\fieldtrip-aleslab-fork
addpath C:\Users\Marlene\Documents\git\ssvepTesting\svndlCopy
addpath C:\Users\Marlene\Documents\git\ssvepTesting\biosemiUpdated
addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFunctions
ft_defaults

% see pred.labels for details
% first 12 are 1st exp, from 12 is the 2nd exp
clearvars;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

% dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};
% dataPath = {'/Users/Marlene/Documents/git/dataLR/LRlongDC/', '/Users/Marlene/Documents/git/dataLR/LRshortDC/'};
dataPath = {'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E2LRlongDC\',...
    'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E1LRshortDC\'};

%%%%%% compute average
for ee=1:2 % which experiment
    clear sbjprediction avPredictions sbj
    load([dataPath{ee} 'predAvRef.mat'])
    load(['avRef-subPredFulE' num2str(ee) '.mat']) % this is for plotting the
%     regression space+time
    
    sbj = sbj(:,1:10);
    avPredictions = averageSbj(sbj);
    avPredictions(1).condLabel = 'originalmotion';
    avPredictions(2).condLabel = 'linear';
    avPredictions(3).condLabel = 'spatial';
    avPredictions(4).condLabel = 'temp';
    avPredictions(5).condLabel = 'spatiotemp';
    avPredictions(6).condLabel = 'SR originalmotion';
    avPredictions(7).condLabel = 'SR linearPred';
    avPredictions(8).condLabel = 'SR spatialPred';
    avPredictions(9).condLabel = 'SR tempPred';
    avPredictions(10).condLabel = 'SR spatiotempPred';
    
    pickElec = [23 10 39]; % Oz PO7 PO8
    nameElec = {'Oz' 'PO7' 'PO8'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% plot average predictions
    figure('Renderer', 'painters', 'Position', [10 10 1300 600])
    condRange = {'LR','SR'};
    for ss=1:2
        for chan=1:length(pickElec)
            subplot(2,3,chan+3*(ss-1)); hold on;
            for mm=2:5 % long/short range
                plot(avPredictions(mm+5*(ss-1)).time,avPredictions(mm+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',1);
            end
            plot(avPredictions(mm+5*(ss-1)).time, squeeze(mean(subPredFul(:,ss,pickElec(chan),:))),'LineWidth',1);
            plot(avPredictions(1+5*(ss-1)).time,avPredictions(1+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',2); % plot AM at the end 
            title([condRange{ss} ' ' nameElec{chan}])
            line([0 400],[0 0],'Color','k','LineStyle','--')
            legend('linear','spatial','temp','S+T','AM','Location','Best')
            if ee==2
                ylim([-4 6])
            elseif ee==1
                ylim([-1.5 2])
            end
        end
    end
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'predictions'],'png')
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'predictions'],'pdf')
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'predictions'],'fig')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% plot interactions
    figure('Renderer', 'painters', 'Position', [10 10 1200 300])
    for chan=1:length(pickElec)
        subplot(1,3,chan); hold on;
        for mm=3:4 % long/short range
            plot(avPredictions(mm).time,avPredictions(2).filteredWave(pickElec(chan),:) - avPredictions(mm).filteredWave(pickElec(chan),:),'LineWidth',2);
            plot(avPredictions(mm+5).time,avPredictions(2+5).filteredWave(pickElec(chan),:) - avPredictions(mm+5).filteredWave(pickElec(chan),:),'LineWidth',2);
        end
        title(['interaction chan' nameElec{chan}])
        line([0 400],[0 0],'Color','k','LineStyle','--')
        legend('spatL','spatS','tempL','tempS','Location','Best')
        if ee==2
            ylim([-2 3])
        elseif ee==1
            ylim([-2 3])
        end
    end
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'interactions'],'png')
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'interactions'],'pdf')
    saveas(gcf,['figures' filesep 'AvRef-E' num2str(ee) 'interactions'],'fig')
    
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
    % RMSE should be difference between prediction - AM but here since it
    % is squared the order is not too important
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
    load(['avRef-rmseCoefE' num2str(ee) '.mat']) % rmseS rmseT rmseS+T x LR/SR 
    load(['avRef-regCoefE' num2str(ee)]); % get the regression coef to get the amount of noise included
    for ss=1:length(sbj)
        nrmseCoef(ss,1) = rmseCoef(ss,3) / (sqrt(rmsNoise(ss,1)^2 + coefF(ss,1,1)*rmsNoise(ss,3)^2 + coefF(ss,1,2)*rmsNoise(ss,4)^2)); % CHECK it is coef*(N^2) not (coef*N)^2
        nrmseCoef(ss,2) = rmseCoef(ss,6) / (sqrt(rmsNoise(ss,6)^2 + coefF(ss,2,1)*rmsNoise(ss,8)^2 + coefF(ss,2,2)*rmsNoise(ss,9)^2 )); 
    end
    
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
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'AvRef-NMRSEboxplot.png'])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'AvRef-NMRSEboxplot.pdf'])
    save(['AvRef-statNRMSE' num2str(ee) '.mat'],'yval') % save NRMSE for stats

end

