% regressions per subj
% same coef for all elec and time points
% linear set to 1

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

% too many variables to clear!!!!
% for ee=1:2 % which experiment

ee=2;

load([dataPath{ee} 'sbjprediction.mat'])
clear percentErrSp percentErrLin percentErrTe percentErrFul subAM subPredFul coefFull correctedSpat correctedTemp
clear rmsTe rmsSp rmsFu subAM correctedSpat correctedTemp subPredFul coefF coefS coefT

for ss=1:length(sbj)
    clear spatInt tempInt subSig regCoefSpa regCoefTem regCoefFul regCoefLin
    clear regCoefS regCoefT regCoefF statsSp statsTe statsFu bintSp bintSp bintFu rSp rTe rFu rintSp rintTe rintFu
    clear coefSp coefTe coefFu predSp predTe predFul
    
    %% get the spatial and temporal interaction terms
    % so substract the linear component from the prediction
    % also get the AM-linear that will be used in the regression
    for numCond=1:2 % long/short
        spatInt(:,:,numCond) = sbj(ss,3+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
        tempInt(:,:,numCond) = sbj(ss,4+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
        subSig(:,:,numCond) = sbj(ss,1+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
    end
    
    % we want to set linear to 1 and only vary the coef for spat and
    % temp components
    % these regressions are therefore done on the residuals of the
    % difference AM-linear (the r2 are for the residuals as well)
    newDim = size(subSig,1)*size(subSig,2);
    for numCond=1:2 % long/short
        regCoefS = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1])]; % as 1*linear + spatial
        regCoefT = [ones(newDim,1) reshape(tempInt(:,:,numCond),[newDim 1])];
        regCoefF = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1]) reshape(tempInt(:,:,numCond),[newDim 1])];
        [coefSp(numCond,:), bintSp(numCond,:,:), rSp(numCond,:), rintSp(numCond,:,:), statsSp(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefS);
        [coefTe(numCond,:), bintSp(numCond,:,:), rTe(numCond,:), rintTe(numCond,:,:), statsTe(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefT);
        [coefFu(numCond,:), bintFu(numCond,:,:), rFu(numCond,:), rintFu(numCond,:,:), statsFu(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefF);
    end
    
    % vector STATS containing, in the following order, the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
    
    % we could calculate r2 from the ones above but might get troubles
    % with the math: can get r2>1
    % instead calculate % error = SSE / sumsqr AM
    % (SSE : sum of squared errors)
    % for linear
    for numCond=1:2 % long/short
        sseLin(ss,numCond,:,:) = sumsqr(sbj(ss,1+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave);
        ssNoise(ss,numCond,:,:) = sumsqr(sbj(ss,1+5*(numCond-1)).data.noiseWave);
        percentErrLin(ss,numCond) = sseLin(ss,numCond,:,:) / ssNoise(ss,numCond,:,:);
        rmsLi(ss,numCond) = rms(sbj(ss,2+5*(numCond-1)).data.filteredWave(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:));
        rmsNoise(ss,numCond,:,:) = rms(sbj(ss,1+5*(numCond-1)).data.noiseWave(:)  );
    end
    %         % for other components
    %         % r2=1-SSE/SStotal
    %         % percentError = 1 - r2
    %         for numCond=1:2 % long/short
    %             percentErrSp(ss,1,numCond) = 1 - statsSp(numCond,1);
    %             percentErrTe(ss,numCond) = 1 - statsTe(numCond,1);
    %             percentErrFu(ss,numCond) = 1 - statsFu(numCond,1);
    %         end
    
    % create sig with the coef
    for numCond=1:2 % long/short
        predSp = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefSp(numCond,2)*spatInt(:,:,numCond);
        sseSp = sumsqr(sbj(ss,1+5*(numCond-1)).data.filteredWave - predSp);
        percentErrSp(ss,numCond) = sseSp / ssNoise(ss,numCond,:,:);
        rmsSp(ss,numCond) = rms(predSp(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:) );
        
        predTe = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefTe(numCond,2)*tempInt(:,:,numCond);
        sseTe = sumsqr(sbj(ss,1+5*(numCond-1)).data.filteredWave - predTe);
        percentErrTe(ss,numCond) = sseTe / ssNoise(ss,numCond,:,:);
        rmsTe(ss,numCond) = rms(predTe(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:) );
        
        predFul = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefFu(numCond,2)*spatInt(:,:,numCond) + coefFu(numCond,3)*tempInt(:,:,numCond);
        sseFul = sumsqr(sbj(ss,1+5*(numCond-1)).data.filteredWave - predFul);
        percentErrFul(ss,numCond) = sseFul / ssNoise(ss,numCond,:,:);
        rmsFu(ss,numCond) = rms( predFul(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:));
        
        correctedSpat(ss,numCond,:,:) = coefFu(numCond,2)*spatInt(:,:,numCond);
        correctedTemp(ss,numCond,:,:) = coefFu(numCond,3)*tempInt(:,:,numCond);
        subPredFul(ss,numCond,:,:) = predFul;
        subAM(ss,numCond,:,:) = sbj(ss,1+5*(numCond-1)).data.filteredWave;
    end
    coefF(ss,:,:) = coefFu(:,2:3);
    coefS(ss,:,:) = coefSp(:,2);
    coefT(ss,:,:) = coefTe(:,2);
end

save(['regCoefE' num2str(ee) '.mat'],'coefF','coefS','coefT');

% need to do it separately for the 2 experiments because the time
% dimensions are different
for numCond=1:2
    avPredFul(numCond,:,:) = squeeze(mean(subPredFul(:,numCond,:,:)));
    avAM(numCond,:,:) = squeeze(mean(subAM(:,numCond,:,:)));
end

rr={'LR','SR'};
figure;
for numCond=1:2
    for chan=1:length(pickElec)
        subplot(2,3,chan+3*(numCond-1)); hold on;
        plot(squeeze(avAM(numCond,pickElec(chan),:)))
        plot(squeeze(avPredFul(numCond,pickElec(chan),:)))
        xlim([1 210])
        legend('AM','pred','Location','Best')
        title(['E' num2str(ee) rr{numCond} 'elec' num2str(pickElec(chan))])
    end
end
saveas(gcf,['figures' filesep 'bestPredElec' num2str(pickElec(chan))],'png')

%     avPercentErrFu(ee,:) = mean(percentErrFul);
%     avPercentErrTe(ee,:) = mean(percentErrTe);
%     avPercentErrSp(ee,:) = mean(percentErrSp);
%     avPercentErrLi(ee,:) = mean(percentErrLin);
%     avCoefFull(ee,:,:) = squeeze(mean(coefF));

%     % only average rms (=RMSE)
%     % rms vs mean over subj
%     avRmsLi(ee,:) = mean(rmsLi );
%     avRmsTe(ee,:) = mean(rmsTe );
%     avRmsSp(ee,:) = mean(rmsSp );
%     avRmsFu(ee,:) = mean(rmsFu );
%     avRmsNoise(ee,:) = mean(rmsNoise );
%     stdRmsLi(ee,:) = std(rmsLi );
%     stdRmsTe(ee,:) = std(rmsTe );
%     stdRmsSp(ee,:) = std(rmsSp );
%     stdRmsFu(ee,:) = std(rmsFu );
%     stdRmsNoise(ee,:) = std(rmsNoise );
rmseCoef = [rmsSp(:,1) rmsTe(:,1) rmsFu(:,1) rmsSp(:,2) rmsTe(:,2) rmsFu(:,2)];
save(['rmseCoefE' num2str(ee) '.mat'],'rmseCoef')

%     % normalise and average rms (=NRMSE)
%     avRmsNormLi(ee,:) = rms(rmsLi ./ rmsNoise);
%     avRmsNormTe(ee,:) = rms(rmsTe ./ rmsNoise);
%     avRmsNormSp(ee,:) = rms(rmsSp ./ rmsNoise);
%     avRmsNormFu(ee,:) = rms(rmsFu ./ rmsNoise);

% looking at average coef is meaningless: depends on the amplitude of the interaction!!!
%%%% plot the interactions
% avCoefFull(ee,numCond,factor)
figure;
for numCond=1:2
    for chan=1:length(pickElec)
        subplot(2,3,chan+3*(numCond-1)); hold on;
        plot(squeeze(mean(correctedSpat(:,numCond,pickElec(chan),:))))
        plot(squeeze(mean(correctedTemp(:,numCond,pickElec(chan),:))))
        ylim([-1.5 1.5])
        title(pickElec(chan))
        xlim([1 210])
        legend('spat','temp')
    end
end
saveas(gcf,['figures' filesep 'fullCoefInteractionsE' num2str(ee)],'png')
% end




% %%%% plot the percentage error
% for ee=1:2
%     figure;
%     bar([avPercentErrLi(ee,1) avPercentErrSp(ee,1) avPercentErrTe(ee,1) avPercentErrFu(ee,1);...
%         avPercentErrLi(ee,2) avPercentErrSp(ee,2) avPercentErrTe(ee,2) avPercentErrFu(ee,2)]);
%     legend('lin','spa','temp','s+t')
%     ylabel('noise error ratio')
%     xticklabels({'LR','SR'})
%     title(['E' num2str(ee) 'with coef'])
%     ylim([0 15])
%     saveas(gcf,['figures' filesep 'noiseErrRatioE' num2str(ee)],'png')
% end
%
% %%%% plot rms...
% for ee=1:2
%     figure;
%     bar([avRmsNormLi(ee,1) avRmsNormSp(ee,1) avRmsNormTe(ee,1) avRmsNormFu(ee,1);...
%         avRmsNormLi(ee,2) avRmsNormSp(ee,2) avRmsNormTe(ee,2) avRmsNormFu(ee,2)]);
%     legend('lin','spa','temp','s+t','Location','Best')
%     ylabel('norm rms')
%     xticklabels({'LR','SR'})
%     title(['E' num2str(ee) 'with coef'])
%     saveas(gcf,['figures' filesep 'coefRmsNormE' num2str(ee)],'png')
% end
% for ee=1:2
%     figure;
%     bar([avRmsNoise(ee,1) avRmsLi(ee,1) avRmsSp(ee,1) avRmsTe(ee,1) avRmsFu(ee,1);...
%         avRmsNoise(ee,2) avRmsLi(ee,2) avRmsSp(ee,2) avRmsTe(ee,2) avRmsFu(ee,2)]);
%     legend('noise','lin','spa','temp','s+t','Location','Best')
%     ylabel('rms')
%     xticklabels({'LR','SR'})
%     title(['E' num2str(ee) 'with coef'])
%     saveas(gcf,['figures' filesep 'coefRmsE' num2str(ee)],'png')
% end


% locX = [1:4 6:9];
% for ee=1:2
%     figure; hold on;
%
%     bar(locX, [avRmsLi(ee,1) avRmsSp(ee,1) avRmsTe(ee,1) avRmsFu(ee,1)...
%         avRmsLi(ee,2) avRmsSp(ee,2) avRmsTe(ee,2) avRmsFu(ee,2)])
%
%     errorbar(locX, [avRmsLi(ee,1) avRmsSp(ee,1) avRmsTe(ee,1) avRmsFu(ee,1)...
%         avRmsLi(ee,2) avRmsSp(ee,2) avRmsTe(ee,2) avRmsFu(ee,2)],...
%     [stdRmsLi(ee,1) stdRmsSp(ee,1) stdRmsTe(ee,1) stdRmsFu(ee,1)...
%         stdRmsLi(ee,2) stdRmsSp(ee,2) stdRmsTe(ee,2) stdRmsFu(ee,2)],'k','Linestyle','none')
% %     legend('noise','lin','spa','temp','s+t','Location','Best')
%     ylabel('RMSE')
% %     xticklabels({'LR','SR'})
%     title(['E' num2str(ee) 'with coef'])
% end
% for ee=1:2
% avRMSE(ee,:) = [avRmsSp(ee,1) avRmsTe(ee,1) avRmsFu(ee,1) avRmsSp(ee,2) avRmsTe(ee,2) avRmsFu(ee,2)];
% stdRMSE(ee,:) = [stdRmsSp(ee,1) stdRmsTe(ee,1) stdRmsFu(ee,1) stdRmsSp(ee,2) stdRmsTe(ee,2) stdRmsFu(ee,2)];
% end


% %%% comparison waveforms AM with Best prediction
% rr={'LR','SR'};
% for chan=1:length(pickElec)
%     figure;
%     for numCond=1:2
%         subplot(2,2,numCond); hold on;
%         plot(squeeze(avAM1(numCond,pickElec(chan),:)))
%         plot(squeeze(avPredFul1(numCond,pickElec(chan),:)))
%         xlim([1 210])
%         legend('AM','pred')
%         title(['EE1' rr{numCond} 'elec' num2str(pickElec(chan))])
%     end
%     for numCond=1:2
%         subplot(2,2,2+numCond); hold on;
%         plot(squeeze(avAM2(numCond,pickElec(chan),:)))
%         plot(squeeze(avPredFul2(numCond,pickElec(chan),:)))
%         xlim([1 210])
%         legend('AM','bestPred')
%         title(['EE2' rr{numCond} 'elec' num2str(pickElec(chan))])
%     end
%     %     saveas(gcf,['figures' filesep 'bestPredElec' num2str(pickElec(chan))],'png')
% end



% 
% %%%% Left over after best fit?
% for numCond=1:2
%     leftOver(numCond,:,:) = squeeze(avAM1(numCond,:,:) - avPredFul1(numCond,:,:));
%     leftOver2(numCond,:,:) = squeeze(avAM2(numCond,:,:) - avPredFul2(numCond,:,:));
% end
% figure;
% for chan=1:length(pickElec)
%     subplot(2,3,chan); hold on;
%     for numCond=1:2
%         plot(squeeze(leftOver(numCond,pickElec(chan),:)))
%         xlim([1 210])
%     end
%     legend('LR','SR')
%     title(['E1elec' num2str(pickElec(chan))])
% end
% for chan=1:length(pickElec)
%     subplot(2,3,chan+3); hold on;
%     for numCond=1:2
%         plot(squeeze(leftOver2(numCond,pickElec(chan),:)))
%         xlim([1 210])
%     end
%     legend('LR','SR')
%     title(['E2elec' num2str(pickElec(chan))])
% end
% figure;plotTopo(squeeze(mean(abs(leftOver(1,:,:)),3)),cfg.layout); colorbar;
% figure;plotTopo(squeeze(mean(abs(leftOver2(1,:,:)),3)),cfg.layout); colorbar;
% figure;plotTopo(squeeze(mean(abs(leftOver(2,:,:)),3)),cfg.layout); colorbar;
% figure;plotTopo(squeeze(mean(abs(leftOver2(2,:,:)),3)),cfg.layout); colorbar;
