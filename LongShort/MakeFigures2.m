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
elec = [23 9 38];

for ee=1:2 % which experiment
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
        for chan=1:length(elec)
            subplot(2,3,chan+3*(ss-1)); hold on;
            for mm=1:5 % long/short range
                plot(avPredictions(mm+5*(ss-1)).time,avPredictions(mm+5*(ss-1)).filteredWave(elec(chan),:),'LineWidth',2);
            end
            legend('AM','linear','spatial','temp','ST','Location','Best')
            title([condRange{ss} num2str(elec(chan))])
            line([0 400],[0 0],'Color','k','LineStyle','--')
        end
    end 
    saveas(gcf,['E' num2str(ee) 'predictions'],'png')
    
    
   
    %%%%%%% RMS root mean square
    % calculate per participant using filteredwave
    for subj=1:length(sbj)
        for chan=1:sbj(1,1).data.nchan
            for ss=1:2 % long/short
                for fact =1:4
                    diffFilt(fact+4*(ss-1),subj,chan,:) = (sbj(subj,1+5*(ss-1)).data.filteredWave(chan,:) - sbj(subj,fact+1+5*(ss-1)).data.filteredWave(chan,:));
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
    
    
    %%% rms bar plot for each sbj
    for chan=1:length(elec)
        figure; hold on;
        for fact=1:8
            bar(fact,mean(rmsFilt(fact,:,elec(chan))))
            scatter(repmat(fact,size(rmsFilt,2),1),rmsFilt(fact,:,elec(chan)))
        end
    end


    %%%%%%% REGRESS for the average signal
    for numCond=1:2 % long/short 
        spatInt(:,:,numCond) = avPredictions(3+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave  ;
        tempInt(:,:,numCond) = avPredictions(4+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave  ;        
        ySig(:,:,numCond) = avPredictions(1+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave ;
    end
    
    % Same coef for time and electrodes
    newDim = [size(ySig,1)*size(ySig,2)];
    for numCond=1:2 % long/short 
        regCoefS = [ones(newDim) reshape(spatInt(:,:,numCond),[newDim 1])]; % as 1*linear + spatial 
        regCoefT = [ones(newDim) reshape(tempInt(:,:,numCond),[newDim 1])]; 
        regCoefF = [ones(newDim) reshape(spatInt(:,:,numCond),[newDim 1]) reshape(tempInt(:,:,numCond),[newDim 1])];
        [coefSp(numCond,:), bintSp(numCond,:), rSp(numCond,:), rintSp(numCond,:), statsSp(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefS);
        [coefTe(numCond,:), bintTe(numCond,:), rTe(numCond,:), rintTe(numCond,:), statsTe(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefT);
        [coefFu(numCond,:), bintFu(numCond,:), rFu(numCond,:), rintFu(numCond,:), statsFu(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefF);
    end
   
    % different coef for different electrodes
    for numCond=1:2 % long/short 
        for chan=2:size(avPredictions(1).filteredWave,1) % first electrode is the reference 0
            regCoefS = [ones(length(ySig),1) spatInt(chan,:,numCond)']; % as 1*linear + spatial 
            regCoefT = [ones(length(ySig),1) tempInt(chan,:,numCond)']; % as 1*linear + spatial 
            regCoefF = [ones(length(ySig),1) spatInt(chan,:,numCond)' tempInt(chan,:,numCond)']; % as 1*linear + spatial 
            [coefSpatial(numCond,chan,:), bintS(numCond,chan,:,:), rS(numCond,chan,:), ...
                rintS(numCond,chan,:,:), statsS(numCond,chan,:)] = ...
                regress(ySig(chan,:,numCond)',regCoefS);
            [coefTemp(numCond,chan,:), bintT(numCond,chan,:,:), rT(numCond,chan,:), ...
                rintT(numCond,chan,:,:), statsT(numCond,chan,:)] = ...
                regress(ySig(chan,:,numCond)',regCoefT);
            [coefFull(numCond,chan,:), bintF(numCond,chan,:,:), rF(numCond,chan,:), ...
                rintF(numCond,chan,:,:), statsF(numCond,chan,:)] = ...
                regress(ySig(chan,:,numCond)',regCoefF);
        end
    end
    % vector STATS containing, in the following order, the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
    
    % topo of r-squares
    figure;
    for numCond=1:2 % long/short
        subplot(3,2,numCond)
        plotTopo(statsS(numCond,:,1),cfg.layout)
        caxis([0 1])
    end
    for numCond=1:2 % long/short
        subplot(3,2,numCond+2)
        plotTopo(statsT(numCond,:,1),cfg.layout)
        caxis([0 1])
    end
    for numCond=1:2 % long/short
        subplot(3,2,numCond+4)
        plotTopo(statsF(numCond,:,1),cfg.layout)
        caxis([0 1])
    end

end




% % reconstruct a "corrected" signal
% for numCond=1:2 % long/short
%     for elec=2:size(avPredictions(1).filteredWave,1) % first electrode is the reference 0
%         bestFit(numCond,elec,:) = coefFull(numCond,elec,2)*avPredictions(2+5*(numCond-1)).filteredWave(elec,:)' + ...
%             coefFull(numCond,elec,3)*spatInt(elec,:,numCond)' + ...
%             coefFull(numCond,elec,4)*tempInt(elec,:,numCond)' ;
%     end
% end
% figure; plot(avPredictions(1).time, avPredictions(1).filteredWave(elec,:));
% hold on; plot(avPredictions(1).time,squeeze(bestFit(1,elec,:)));
% figure; plot(avPredictions(6).time, avPredictions(6).filteredWave(elec,:));
% hold on; plot(avPredictions(6).time,squeeze(bestFit(2,elec,:)));
% 
% % left over after best-fitting
% for numCond=1:2 % long/short
%     leftOver(numCond,:,:) = avPredictions(1+5*(numCond-1)).filteredWave - squeeze(bestFit(numCond,:,:)) ;
% end
% figure; plot(avPredictions(1).time, avPredictions(1).filteredWave(elec,:)-squeeze(bestFit(1,elec,:))');
% figure; plot(avPredictions(1).time, avPredictions(6).filteredWave(elec,:)-squeeze(bestFit(2,elec,:))');
% figure;
% for numCond=1:2 % long/short
%     subplot(1,2,numCond)
%     plotTopo(squeeze(mean(abs(leftOver(numCond,:,:)),3)),cfg.layout)
%     colorbar;
% end
    