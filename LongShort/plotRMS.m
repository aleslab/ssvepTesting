
clearvars;
close all;

for ee=1:2
    load(['allRMSe' num2str(ee)]);
    load(['regCoefE' num2str(ee)]);
    figure; hold on
    subplot(2,2,1);hold on;
    bar([mean(rmsSignal(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))])
    errorbar([mean(rmsSignal(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))], ...
        [std(rmsSignal(:,1:2)) std(rmseNoCoef(:,1:4)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['rmsE' num2str(ee) 'LR'])
    subplot(2,2,2);hold on;
    bar([mean(rmsSignal(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))])
    errorbar([mean(rmsSignal(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))], ...
        [std(rmsSignal(:,6:7)) std(rmseNoCoef(:,5:8)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
    title(['rmsE' num2str(ee) 'SR'])
    % x-axis: AM lin rmseNoCoef (lin s t s+t) rmseCoef (s t s+t)
    
    subplot(2,2,3);hold on;
    bar(mean(rmsNoise))
    errorbar(mean(rmsNoise), std(rmsNoise),'k','Linestyle','none')
    title(['raw noiseE' num2str(ee)])
    ylabel('rms noise per cond')
    % x-axis: LR for AM lin S T ST, SR for AM lin S T ST
    
    % correct noise
    % order: linear spatial temporal s+t
    % eg noiseLin = sqrt(noiseAM*noiseAM + noiseLin*noiseLin)
    % should be similar to noiseAM * sqrt(3)
    % the std per participant for the noise is around 0.5 calculated by:
    % sqrt([(noisecorrected / sqrt(2/pi))^2] * (4-pi)/2)
    noiseCorr = zeros(length(rmsNoise),11);
    noiseCorr2 = zeros(length(rmsNoise),11); % this is just to have a look with the sqrt results
    %     for ss=1:length(rmsNoise)
    %         for motrange=1:2
    %             noiseCorr(ss,1+4*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2); % AM + lin
    %             noiseCorr2(ss,1+4*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3); % AM + lin
    %             noiseCorr(ss,2+4*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2); % AM + lin + spat
    %             noiseCorr2(ss,2+4*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(4);
    %             noiseCorr(ss,3+4*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % AM + lin + temp
    %             noiseCorr2(ss,3+4*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(4);
    %             noiseCorr(ss,4+4*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % AM + lin + spat + temp
    %             noiseCorr2(ss,4+4*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(5);
    %        end
    %     end
    %%% not including AM noise for the rmse
    for ss=1:length(rmsNoise)
        for motrange=1:2
            noiseCorr(ss,1+7*(motrange-1)) = sqrt(rmsNoise(ss,1+5*(motrange-1))^2 + rmsNoise(ss,2+5*(motrange-1))^2); % AM + lin
            noiseCorr2(ss,1+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3); % AM + lin
            noiseCorr(ss,2+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2); % lin + spat
            noiseCorr2(ss,2+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,3+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % lin + temp
            noiseCorr2(ss,3+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(3);
            noiseCorr(ss,4+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + rmsNoise(ss,3+5*(motrange-1))^2 + rmsNoise(ss,4+5*(motrange-1))^2); % lin + spat + temp
            noiseCorr2(ss,4+7*(motrange-1)) = rmsNoise(ss,1+5*(motrange-1)) * sqrt(4);
            % noise from regression uses the coef
            noiseCorr(ss,5+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefS(ss,motrange)*rmsNoise(ss,3+5*(motrange-1))^2); 
            noiseCorr(ss,6+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefT(ss,motrange)*rmsNoise(ss,4+5*(motrange-1))^2); 
            noiseCorr(ss,7+7*(motrange-1)) = sqrt(rmsNoise(ss,2+5*(motrange-1))^2 + coefF(ss,motrange,1)*rmsNoise(ss,3+5*(motrange-1))^2 + coefF(ss,motrange,2)*rmsNoise(ss,4+5*(motrange-1))^2);             
        end
    end
    subplot(2,2,4);hold on;
    bar(mean(noiseCorr))
    errorbar(mean(noiseCorr), std(noiseCorr),'k','Linestyle','none')
    title(['corrected noise no coef E' num2str(ee)])
    ylabel('rms noise')
    
    % apply normalisation by noise for each participant
%     normRMSlr = [rmsSignal(:,1:2)./rmsNoise(:,1:2) rmseNoCoef(:,1:4) ./ noiseCorr(:,1:4)];
%     normRMSsr = [rmsSignal(:,6:7)./rmsNoise(:,6:7) rmseNoCoef(:,5:8) ./ noiseCorr(:,5:8)];
    normRMSlr = [rmsSignal(:,1:2)./rmsNoise(:,1:2) rmseNoCoef(:,1:4) ./ noiseCorr(:,1:4) rmseCoef(:,1:3) ./ noiseCorr(:,5:7)];
    normRMSsr = [rmsSignal(:,6:7)./rmsNoise(:,6:7) rmseNoCoef(:,5:8) ./ noiseCorr(:,8:11) rmseCoef(:,4:6) ./ noiseCorr(:,12:14)];
    
    figure; hold on;
    subplot(2,1,1);hold on;
    bar(mean(normRMSlr))
    errorbar(mean(normRMSlr),std(normRMSlr),'k','Linestyle','none')
    ylabel('RMS signal/ rms noise')
    line([1 9],[1 1],'Color','r','LineWidth',3)
    title(['E' num2str(ee) ' LR'])
    subplot(2,1,2);hold on;
    bar(mean(normRMSsr))
    errorbar(mean(normRMSsr),std(normRMSsr),'k','Linestyle','none')
    ylabel('RMS signal/ rms noise')
    line([1 9],[1 1],'Color','r','LineWidth',3)
    title(['E' num2str(ee) ' SR'])
    saveas(gcf,['figures' filesep 'RMScorrected' num2str(ee)],'png')
    
    %     norm2RMSlr = [rmsSignal(:,1:2)./rmsNoise(:,1:2) rmseNoCoef(:,1:4) ./ noiseCorr2(:,1:4)];
    %     norm2RMSsr = [rmsSignal(:,6:7)./rmsNoise(:,6:7) rmseNoCoef(:,5:8) ./ noiseCorr2(:,5:8)];
    %
    %     figure; hold on;
    %     subplot(2,1,1);hold on;
    %     bar(mean(norm2RMSlr))
    %     errorbar(mean(norm2RMSlr),std(norm2RMSlr),'k','Linestyle','none')
    %     ylabel('RMS signal/ rms noise')
    %     line([1 6],[1 1],'Color','r','LineWidth',3)
    %     title(['E' num2str(ee) ' LR'])
    %     subplot(2,1,2);hold on;
    %     bar(mean(norm2RMSsr))
    %     errorbar(mean(norm2RMSsr),std(norm2RMSsr),'k','Linestyle','none')
    %     ylabel('RMS signal/ rms noise')
    %     line([1 6],[1 1],'Color','r','LineWidth',3)
    %     title(['E' num2str(ee) ' SR'])
    
end


