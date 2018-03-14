

% create a new field "filteredWave" after low-pass filter
% warning about filter is NORMAL
for condIdx=1:length(gpCorrection)
    filtIdx = determineFilterIndices( 'low49', gpCorrection(condIdx).freq, gpCorrection(condIdx).i1f1 );
    
    %Create a logical matrix selecting frequency components.
    filtMat = false(size(gpCorrection(condIdx).amp));
    filtMat(:,filtIdx) = true;
    
    %Combine the filter and sig vaules with a logical AND.
    % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
    gpCorrection(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
    cfg.activeFreq =  gpCorrection(condIdx).activeFreq;
    
    [ gpCorrection(condIdx).filteredWave ] = filterSteadyState( cfg, gpCorrection(condIdx) );
end



% waveforms

chan = [23 19]; 
chanName = {'Oz' 'Pz'};

limERPmin = -1.2;
limERPmax = 1.2;

colors = {'r','k','c','g','b'};

% because of the pb with pre-processing starting at 1.529, the stimulus
% appears on the left only 80 ms later. Correction = 1st time point is 42
% gpCorrection(condIdx).time([42:204 1:41])

for plotCh=1:length(chan)
channel = chan(plotCh);
figure; hold on;
for condIdx=1:5
    plot(gpCorrection(condIdx).time,gpCorrection(condIdx).filteredWave(channel,[42:204 1:41]),colors{condIdx},'LineWidth',2)
    ylim([limERPmin limERPmax])
end
xlabel('time (ms)')
legend('M','L','L+S','L+T','L+S+T','Location','NorthWest')
title(['LR ' chanName{plotCh}])
saveas(gcf,['LR_' chanName{plotCh} '.jpg'])
saveas(gcf,['LR_' chanName{plotCh} '.epsc'])
saveas(gcf,['LR_' chanName{plotCh} '.m'])

figure; hold on;
for condIdx=6:10
    plot(gpCorrection(condIdx).time,gpCorrection(condIdx).filteredWave(channel,[42:204 1:41]),colors{condIdx-5},'LineWidth',2)
    ylim([limERPmin limERPmax])
end
legend('M','L','L+S','L+T','L+S+T','Location','NorthWest')
xlabel('time (ms)')
title(['SR ' chanName{plotCh}])
saveas(gcf,['SR_' chanName{plotCh} '.jpg'])
saveas(gcf,['SR_' chanName{plotCh} '.epsc'])
saveas(gcf,['SR_' chanName{plotCh} '.m'])

end



% power spectrum
limFq = 0.8;
plotName = {'L','L+S','L+T','L+S+T'};
for plotCh=1:length(chan)
    channel = chan(plotCh);

figure; hold on
for condition=1:4
    subplot(2,2,condition);hold on;
    bar(squeeze(gpCorrection(1).freq((gpCorrection(condition).freq<15.5 & gpCorrection(1).freq>0))),squeeze(gpCorrection(1).amp(channel,2:31)),0.8,'FaceColor','r')
    bar(squeeze(gpCorrection(condition+1).freq((gpCorrection(condition).freq<15.5 & gpCorrection(condition).freq>0))),squeeze(gpCorrection(condition+1).amp(channel,2:31)),0.3,'FaceColor','b')
    ylim([0 limFq])
    title([plotName{condition} ' ' chanName{plotCh}])
end
saveas(gcf,['LRfreq_' chanName{plotCh} '.jpg'])

figure; hold on
for condition=6:9
    subplot(2,2,condition-5);hold on;
    bar(squeeze(gpCorrection(1).freq((gpCorrection(condition).freq<15.5 & gpCorrection(1).freq>0))),squeeze(gpCorrection(1).amp(channel,2:31)),0.8,'FaceColor','r')
    bar(squeeze(gpCorrection(condition+1).freq((gpCorrection(condition).freq<15.5 & gpCorrection(condition).freq>0))),squeeze(gpCorrection(condition+1).amp(channel,2:31)),0.3,'FaceColor','b')
    ylim([0 limFq])
    title([plotName{condition-5} ' ' chanName{plotCh}])
end
saveas(gcf,['SRfreq_' chanName{plotCh} '.jpg'])

end