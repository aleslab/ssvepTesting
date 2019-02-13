% regressions

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
    
    %%%%%%% REGRESS for the average signal
    for numCond=1:2 % long/short 
        spatInt(:,:,numCond) = avPredictions(3+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave  ;
        tempInt(:,:,numCond) = avPredictions(4+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave  ;        
        ySig(:,:,numCond) = avPredictions(1+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave ;
    end
    
    % Same coef for time and electrodes
    newDim = size(ySig,1)*size(ySig,2);
    for numCond=1:2 % long/short 
        regCoefS = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1])]; % as 1*linear + spatial 
        regCoefT = [ones(newDim,1) reshape(tempInt(:,:,numCond),[newDim 1])]; 
        regCoefF = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1]) reshape(tempInt(:,:,numCond),[newDim 1])];
        [coefSp(numCond,:), bintSp(numCond,:,:), rSp(numCond,:), rintSp(numCond,:,:), statsSp(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefS);
        [coefTe(numCond,:), bintTe(numCond,:,:), rTe(numCond,:), rintTe(numCond,:,:), statsTe(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefT);
        [coefFu(numCond,:), bintFu(numCond,:,:), rFu(numCond,:), rintFu(numCond,:,:), statsFu(numCond,:)] = regress(reshape(ySig(:,:,numCond),[newDim 1]),regCoefF);
    end
    
    for numCond=1:2 % long/short
        sseLin = sumsqr(avPredictions(1+5*(numCond-1)).filteredWave - avPredictions(2+5*(numCond-1)).filteredWave);
        ssAM = sumsqr(avPredictions(1+5*(numCond-1)).filteredWave);
        percentErrLin(numCond) = sseLin / ssAM;
    end
    % for other components
    % r2=1-SSE/SStotal
    % percentError = 1 - r2
    for numCond=1:2 % long/short
        percentErrSp(numCond) = 1 - statsSp(numCond,1);
        percentErrTe(numCond) = 1 - statsTe(numCond,1);
        percentErrFu(numCond) = 1 - statsFu(numCond,1);
        rsqSp(numCond) = statsSp(numCond,1);
        rsqTe(numCond) = statsTe(numCond,1);
        rsqFu(numCond) = statsFu(numCond,1);
    end
        
    for chan=1:length(pickElec)
        figure;
        for numCond=1:2 % long/short
            subplot(2,2,numCond+1*(numCond-1)); hold on;
            plot(avPredictions(1+5*(numCond-1)).filteredWave(pickElec(chan),:))
            plot(avPredictions(2+5*(numCond-1)).filteredWave(pickElec(chan),:))
            legend('AM','Linear')
            title([condRange{numCond} num2str(pickElec(chan))])
            subplot(2,2,numCond+1+1*(numCond-1)); hold on;
            plot(coefFu(numCond,2)*spatInt(pickElec(chan),:,numCond))
            plot(coefFu(numCond,3)*tempInt(pickElec(chan),:,numCond))
            legend('Spatial','Temporal')
            title([condRange{numCond} num2str(pickElec(chan))])
        end
    end
   
   
   clear regCoefS regCoefT regCoefF   
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
    