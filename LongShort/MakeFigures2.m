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
    clear diffFilt rmsFilt
    load([dataPath{ee} 'sbjprediction.mat'])

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
            legend('AM','linear','spatial','temp','ST','Location','Best')
            title([condRange{ss} num2str(pickElec(chan))])
            line([0 400],[0 0],'Color','k','LineStyle','--')
        end
    end 
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'predictions'],'png')
    
    
   
    %%%%%%% RMS root mean square
    % calculate per participant using filteredwave
    for subj=1:length(sbj)
        for chan=1:sbj(1,1).data.nchan
            for ss=1:2 % long/short
                for fact =1:4
                    diffFilt(fact+4*(ss-1),subj,chan,:) = (sbj(subj,1+5*(ss-1)).data.filteredWave(chan,:) - sbj(subj,fact+1+5*(ss-1)).data.filteredWave(chan,:));
                    rmsFilt(fact+4*(ss-1),subj,chan) = rms(diffFilt(fact+4*(ss-1),subj,chan,:));
                end
                rmsNoise(ss,subj,chan) = rms(sbj(subj,1+5*(ss-1)).data.noiseWave(chan,:));
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
    
    %%% rms of sbj across electrodes
    rmsNoiseAll = rms(rms(rmsNoise(:,:,:),3),2);
    allRMS = rms(rms(rmsFilt(:,:,:),3),2);
    newRMS = reshape(allRMS,[4 2]);
    figure;
    bar([rmsNoiseAll'; newRMS]')  
    legend('noise','linear','spatial','temporal','s+t','Location','Best')
    xticklabels({'LR','SR'})
    ylabel('rms with rms(sbj)')
    title(['E' num2str(ee)])
    ylim([0 1.4])
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'rms'],'png')
    
end




