

dataIn = '/Users/marleneponcet/Documents/data/seqContext/cleanData/';

% find sequence > 2 for each condition
allCond = unique(cleanData.trialinfo(:,1));
for condition = 1:length(allCond)
    cfg = [];
    cfg.trials = find(cleanData.trialinfo(:,1)==allCond(condition) & cleanData.trialinfo(:,4)>2);
    cfg.preproc.baselinewindow  = [-0.2 0];
    cfg.demean = 'no'; % no baseline correction
    avData(condition) = ft_timelockanalysis(cfg, cleanData);
end

figure; hold on
for condition = 1:4
    plot(avData(condition).time,avData(condition).avg(23,:),'LineWidth',2)
end
timePerStim = length(avData(1).time)/6;
stimTime = avData(1).time(1:timePerStim:length(avData(1).time));
for ll=2:length(stimTime)
    line([stimTime(ll) stimTime(ll)],[-6 10],'Color','k','LineStyle','--')
end
legend('motion','seq','reversal','two')


% try to look only at last stim response
for cond=1:4
    cutData(cond,:,:) = avData(cond).avg(:,288:end) - mean(avData(cond).avg(:,144:287),2);
end
figure; hold on
for condition = 1:4
    plot(squeeze(cutData(condition,23,:)),'LineWidth',2)
end
