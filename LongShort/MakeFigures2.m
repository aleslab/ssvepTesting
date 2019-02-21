% Make figures
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
pickElec = [23 9 38];

for ee=1:2 % which experiment
    clear diffFilt rmsFilt rmseNoCoef rmsNoise rmsSignal rmseCoef rmseNoCoef testrmsNoise
    load([dataPath{ee} 'sbjprediction.mat'])
    load(['rmseCoefE' num2str(ee) '.mat'])

    %%%%%% compute average
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
    avPredictions(11).condLabel = 'LR STnl';
    avPredictions(12).condLabel = 'SR STnl';
    
    %%%%%% plot average predictions
    figure('Renderer', 'painters', 'Position', [10 10 1200 500])
    condRange = {'LR','SR'};
    for ss=1:2
        for chan=1:length(pickElec)
            subplot(2,3,chan+3*(ss-1)); hold on;
            for mm=1:5 % long/short range
                plot(avPredictions(mm+5*(ss-1)).time,avPredictions(mm+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',2);
            end
            title([condRange{ss} num2str(pickElec(chan))])
            line([0 400],[0 0],'Color','k','LineStyle','--')
            legend('AM','linear','spatial','temp','ST','Location','Best')
       end
    end 
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'predictions'],'png')
    
    
   
    %%%%%%% RMS root mean square
    % calculate per participant using filteredwave
    for subj=1:length(sbj)
        for chan=1:sbj(1,1).data.nchan
            for ss=1:2 % long/short
                for fact =1:4
                    diffFilt(fact+4*(ss-1),subj,chan,:) = (sbj(subj,fact+1+5*(ss-1)).data.filteredWave(chan,:) - sbj(subj,1+5*(ss-1)).data.filteredWave(chan,:));
                    rmsFilt(fact+4*(ss-1),subj,chan) = rms(diffFilt(fact+4*(ss-1),subj,chan,:));
                end
            end
        end
    end
    
    %%% rms topography on averaged signal
    figure;hold on
    for fact=1:8
        tmp = squeeze(mean(diffFilt(fact,:,:,:),2));
        subplot(2,4,fact)
        plotTopo(rms(tmp,2),cfg.layout)
        colorbar
    end
    
    
%     %%% rms bar plot for each sbj
%     for chan=1:length(pickElec)
%         figure; hold on;
%         for fact=1:8
%             bar(fact,rms(rmsFilt(fact,:,pickElec(chan))))
%             scatter(repmat(fact,size(rmsFilt,2),1),rmsFilt(fact,:,pickElec(chan)))
%         end
%     end
    
%     %%% rms of sbj across electrodes
%     rmsNoiseAll = rms(rms(rmsNoise(:,:,:),3),2);
%     stdNoiseAll = std(rms(rmsNoise(:,:,:),3),0,2);
%     allRMS = rms(rms(rmsFilt(:,:,:),3),2);
%     stdAllRMS = std(rms(rmsFilt(:,:,:),3),0,2);
%     
%     newRMS = reshape(allRMS,[4 2]); newSTD = reshape(stdAllRMS,[4 2]);
%     newRMScat = [rmsNoiseAll'; newRMS];  newSTDcat = [stdNoiseAll'; newSTD];
%     newRMScat = reshape(newRMScat,10,1); newSTDcat = reshape(newSTDcat,10,1);
%     
%     colBar = {'b','m','c','g','r','b','m','c','g','r'};
%     locX = [1:5 7:11];
%     figure; hold on;
%     for dd=1:10
%         bar(locX(dd),newRMScat(dd),colBar{dd})
%     end
%     errorbar(locX, newRMScat, newSTDcat,'k','Linestyle','none')
%     legend('noise','linear','spatial','temporal','s+t','Location','Best')
%     xticklabels({'','LR','','','SR'})
%     ylabel('rms with rms(sbj)')
%     title(['E' num2str(ee)])
%     ylim([0 2])
%     saveas(gcf,['figures' filesep 'E' num2str(ee) 'rms'],'png')

    %%%%%%%%%%%%%%%%%
    % look at rms per condition
    for ss=1:length(sbj)
        for numCond=1:10
            testrmsNoise(ss,numCond,:) = rms(sbj(ss,numCond).data.noiseWave(:)  );
            rmsSignal(ss,numCond,:) = rms(sbj(ss,numCond).data.filteredWave(:)  );
        end
    end
    figure;hold on;
    bar(mean(testrmsNoise))
    errorbar(mean(testrmsNoise), std(testrmsNoise),'k','Linestyle','none')
    title(['E' num2str(ee)])
    ylabel('rms noise per cond')
    saveas(gcf,['figures' filesep 'RMSnoiseAllCond' num2str(ee)],'png')
   
    figure;hold on;
    bar(mean(rmsSignal))
    errorbar(mean(rmsSignal), std(rmsSignal),'k','Linestyle','none')
    title(['E' num2str(ee)])
    ylabel('rms Signal per cond')
    
    %%%%%%%
    % get all the rms we need across electrodes
    rmseNoCoef(:,:) = rms(rmsFilt(:,:,:),3)';
    rmsNoise = rms(testrmsNoise(:,:,:),3);
    save(['allRMSe' num2str(ee) '.mat'],'rmsSignal','rmseNoCoef','rmsNoise','rmseCoef')
%     rmseAvg(ee,:) = allRMS';
%     rmseStd(ee,:) = stdAllRMS';
%     rmseAvgNoise(ee,:) = mean(testrmsNoise);
%     rmseStdNoise(ee,:) = std(testrmsNoise);    
%     rmsAvg(ee,:) = mean(rmsSignal);
%     rmsStd(ee,:) = std(rmsSignal);        
end





figure;
for ee=1:2
    load(['allRMSe' num2str(ee)]);
    subplot(2,2,1+2*(ee-1));hold on;
    bar([mean(rmsSignal(:,1:2)) mean(rmsNoise(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))])
    errorbar([mean(rmsSignal(:,1:2)) mean(rmsNoise(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))], ...
        [std(rmsSignal(:,1:2)) std(rmsNoise(:,1:2)) std(rmseNoCoef(:,1:4)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['E' num2str(ee) 'LR'])
    subplot(2,2,2+2*(ee-1));hold on;
    bar([mean(rmsSignal(:,6:7)) mean(rmsNoise(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))])
    errorbar([mean(rmsSignal(:,6:7)) mean(rmsNoise(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))], ...
        [std(rmsSignal(:,6:7)) std(rmsNoise(:,6:7)) std(rmseNoCoef(:,5:8)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['E' num2str(ee) 'SR'])
end
saveas(gcf,['figures' filesep 'RMSallCond' num2str(ee)],'png')




