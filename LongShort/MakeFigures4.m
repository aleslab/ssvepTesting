% Make figures
% predictions (with AM)
% Interactions
% (for RMS see RMSanalysis)

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
dataPath = {'/Users/Marlene/Documents/git/dataLR/LRlongDC/', '/Users/Marlene/Documents/git/dataLR/LRshortDC/'};


%%%%%% compute average
for ee=1:2 % which experiment
    clear sbjprediction avPredictions
    load([dataPath{ee} 'sbjprediction.mat'])
    load(['subPredFulE' num2str(ee) '.mat'])
    
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
            legend('linear','spatial','temp','S+T','regress','AM','Location','Best')
            if ee==2
                ylim([-4 6])
            elseif ee==1
                ylim([-1.5 2])
            end
        end
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'predictions'],'png')
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'predictions'],'fig')
    
    
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
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'interactions'],'png')
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'interactions'],'fig')
    
end

