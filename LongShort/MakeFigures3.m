% Make figures
% topo of SSVEP
% predictions (with AM)
% Interactions
% (for RMS see plotRMS)

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


%%%%%% compute average
for ee=1:2 % which experiment
    clear sbjprediction avPredictions
    load([dataPath{ee} 'sbjprediction.mat'])
    load(['rmseCoefE' num2str(ee) '.mat'])
    
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Topographies to determine which electrodes to include
    condition = [1 6]; condName = {'long range' 'short range'};
    harm = avPredictions(1).freq(avPredictions(1).i1f1).*[1:6];
    for hh=1:length(harm)
        harmIndex(hh) = find(avPredictions(1).freq==harm(hh));
    end
    
    % for the short and long range conditions = noisy
    figure('Renderer', 'painters', 'Position', [10 10 1200 800])
    for cc=1:length(condition)
        for hh=1:4 % length(harmIndex)
            subplot(2,4,hh+4*(cc-1)); hold on;
            plotTopo(avPredictions(condition(cc)).amp(:,harmIndex(hh)),cfg.layout);
            maxVal(cc,hh) = max(avPredictions(condition(cc)).amp(:,harmIndex(hh)));
            maxElec(cc,hh) = find(avPredictions(condition(cc)).amp(:,harmIndex(hh)) == maxVal(cc,hh));
            colorbar;
            title([num2str(harm(hh)) 'Hz'])
        end
    end
    colormap(jmaColors('hotcortex'));
    colormap('hot');
    saveas(gcf,['figures' filesep 'topoE' num2str(ee) '.png'])
    
    
    pickElec = [23 126 38];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            if ee==2
                ylim([-4 6])
            elseif ee==1
                ylim([-1.5 2])
            end
        end
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'predictions'],'png')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% plot interactions
    %     figure('Renderer', 'painters', 'Position', [10 10 1200 500])
    %     for ss=1:2
    %         for chan=1:length(pickElec)
    %             subplot(2,3,chan+3*(ss-1)); hold on;
    %             for mm=3:5 % long/short range
    %                 plot(avPredictions(mm+5*(ss-1)).time,avPredictions(2+5*(ss-1)).filteredWave(pickElec(chan),:) - avPredictions(mm+5*(ss-1)).filteredWave(pickElec(chan),:),'LineWidth',2);
    %             end
    %             title(['interaction' condRange{ss} num2str(pickElec(chan))])
    %             line([0 400],[0 0],'Color','k','LineStyle','--')
    %             legend('spatial','temp','ST','Location','Best')
    %             if ee==2
    %                 ylim([-4 6])
    %             elseif ee==1
    %                 ylim([-2 3])
    %             end
    %        end
    %     end
    %     saveas(gcf,['figures' filesep 'E' num2str(ee) 'interactions'],'png')
    figure('Renderer', 'painters', 'Position', [10 10 1200 300])
    for chan=1:length(pickElec)
        subplot(1,3,chan); hold on;
        for mm=3:4 % long/short range
            plot(avPredictions(mm).time,avPredictions(2).filteredWave(pickElec(chan),:) - avPredictions(mm).filteredWave(pickElec(chan),:),'LineWidth',2);
            plot(avPredictions(mm+5).time,avPredictions(2+5).filteredWave(pickElec(chan),:) - avPredictions(mm+5).filteredWave(pickElec(chan),:),'LineWidth',2);
        end
        title(['interaction chan' num2str(pickElec(chan))])
        line([0 400],[0 0],'Color','k','LineStyle','--')
        legend('spatL','spatS','tempL','tempS','Location','Best')
        if ee==2
            ylim([-2 3])
        elseif ee==1
            ylim([-2 3])
        end
    end
    saveas(gcf,['figures' filesep 'E' num2str(ee) 'interactions'],'png')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % per sbj
        
        figure('Renderer', 'painters', 'Position', [10 10 1200 300])
        for chan=1:length(pickElec)
            subplot(1,3,chan); hold on;
            plot(avPredictions(1).time,avPredictions(1).filteredWave(pickElec(chan),:)-avPredictions(4).filteredWave(pickElec(chan),:),'LineWidth',2)
            plot(avPredictions(1).time,avPredictions(6).filteredWave(pickElec(chan),:)-avPredictions(9).filteredWave(pickElec(chan),:),'LineWidth',2)
            line([0 400],[0 0],'Color','k','LineStyle','--')
            legend('LR','SR')
            title(['unexplained chan' num2str(pickElec(chan))])
            ylim([-2 3])
        end
        saveas(gcf,['figures' filesep 'Rest3chan E' num2str(ee)],'png')
        
%         figure;
%         for ss=1:length(sbj)
%             subplot(4,4,ss); hold on;
%             plot(sbj(ss,1).data.time,sbj(ss,1).data.filteredWave(pickElec,:)-sbj(ss,4).data.filteredWave(pickElec,:),'LineWidth',2);
%             plot(sbj(ss,1).data.time,sbj(ss,6).data.filteredWave(pickElec,:)-sbj(ss,9).data.filteredWave(pickElec,:),'LineWidth',2);
%         end       
        
%         figure;
%         for ss=1:length(sbj)
%             subplot(4,4,ss); hold on;
%             plot(sbj(ss,1).data.time,sbj(ss,1).data.filteredWave(pickElec,:)-sbj(ss,4).data.filteredWave(pickElec,:),'LineWidth',2);
%     %         plot(avPredictions(1).time,avPredictions(1).filteredWave(pickElec,:) - avPredictions(4).filteredWave(pickElec,:),'LineWidth',2);
%             plot(sbj(ss,1).data.time,sbj(1).data.noiseWave(pickElec,:) ,'LineWidth',2);
%             testRMS(ss) = rms(sbj(ss,1).data.filteredWave(pickElec,:)-sbj(ss,4).data.filteredWave(pickElec,:));
%             title(['rms' num2str(testRMS(ss))])
%         end
%         figure;
%         for ss=1:length(sbj)
%             subplot(4,4,ss); hold on;
%             plot(sbj(ss,1).data.time,sbj(ss,6).data.filteredWave(pickElec,:)-sbj(ss,9).data.filteredWave(pickElec,:),'LineWidth',2);
%     %         plot(avPredictions(1).time,avPredictions(6).filteredWave(pickElec,:) - avPredictions(9).filteredWave(pickElec,:),'LineWidth',2);
%             plot(sbj(ss,6).data.time,sbj(6).data.noiseWave(pickElec,:) ,'LineWidth',2);
%             plot(sbj(ss,6).data.time,rmsTopoNoise(ss,5,9),'LineWidth',2);
%             testRMSShort(ss) = rms(sbj(ss,6).data.filteredWave(pickElec,:)-sbj(ss,9).data.filteredWave(pickElec,:));
%             title(['rms' num2str(testRMSShort(ss))])
%         end
%     
%         pickElec = 16;
%         figure; hold on
%         subplot(2,1,1); hold on;
%         plot(avPredictions(1).time,avPredictions(1).filteredWave(pickElec,:)-avPredictions(2).filteredWave(pickElec,:),'LineWidth',2);
%         plot(avPredictions(1).time,avPredictions(1).filteredWave(pickElec,:)-avPredictions(4).filteredWave(pickElec,:),'LineWidth',2);
%         plot(avPredictions(1).time,avPredictions(1).noiseWave(pickElec,:),'LineWidth',2);
%         line([0 400],[0 0],'Color','k','LineStyle','--')
%         title(['LRe' num2str(ee)])
%         subplot(2,1,2); hold on;
%         plot(avPredictions(1).time,avPredictions(6).filteredWave(pickElec,:)-avPredictions(7).filteredWave(pickElec,:),'LineWidth',2);
%         plot(avPredictions(1).time,avPredictions(6).filteredWave(pickElec,:)-avPredictions(9).filteredWave(pickElec,:),'LineWidth',2);
%         plot(avPredictions(1).time,avPredictions(6).noiseWave(pickElec,:),'LineWidth',2);
%         line([0 400],[0 0],'Color','k','LineStyle','--')
%         legend('AM-lin','AM-temp','noiseAM','location','best')
%         title(['SRe' num2str(ee)])
%         saveas(gcf,['figures' filesep 'chan16 E' num2str(ee)],'png')
    
end

