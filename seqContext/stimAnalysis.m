
clearvars;
% filesDir = 'C:\Users\Marlene\Documents\dataStAndrews\seqContext\clean\';
filesDir = '/Users/marleneponcet/Documents/data/seqContext/clean/';
filesList = dir([filesDir '*stim.mat']);

plotElec = [23 9];

for ss=1:length(filesList)
clear cleanData
load([filesDir filesList(ss).name])
% find sequence > 2 for each condition
allCond = unique(cleanData.trialinfo(:,1));
for condition = 1:length(allCond)
    cfg = [];
    cfg.trials = find(cleanData.trialinfo(:,1)==allCond(condition) & cleanData.trialinfo(:,4)>2);
    cfg.preproc.demean = 'no'; % no baseline correction
    cleanData.elec=elec;
    avData(ss,condition) = ft_timelockanalysis(cfg, cleanData);
end
% for condition = 1:length(allCond)
%    	index = find(cleanData.trialinfo(:,1)==allCond(condition) & cleanData.trialinfo(:,4)>2);
%     clear condTrial
%     for tt=1:length(index)
%         condTrial(tt,:,:) = cleanData.trial{1,index(tt)};
%     end
%     avDataT(:,:,condition) = squeeze(mean(condTrial));
% end

for ee=1:length(plotElec)
figure; hold on
for condition = 1:length(allCond)
    plot(avData(ss,condition).time,avData(ss,condition).avg(plotElec(ee),:),'LineWidth',2)
end
legend('motion','seq','reversal','two','rand')
title(['S' num2str(ss) ' E' num2str(plotElec(ee))])
saveas(gcf,['figures', filesep 'S' num2str(ss) ' E' num2str(plotElec(ee))],'png')
end


% figure; hold on
% for condition = 1:length(allCond)
%     plot(avDataT(23,:,condition),'LineWidth',2)
% end

end




for condition = 1:length(allCond)
    grandavg(condition) = ft_timelockgrandaverage(cfg,avData(1,condition),avData(2,condition),avData(3,condition));
end

for ee=1:length(plotElec)
figure; hold on
for condition = 1:length(allCond)
    plot(grandavg(condition).time,grandavg(condition).avg(plotElec(ee),:),'LineWidth',2)
end
legend('motion','seq','reversal','two','rand')
title(['Avg' ' E' num2str(plotElec(ee))])
saveas(gcf,['figures', filesep 'Avg' ' E' num2str(plotElec(ee))],'png')
end

ft_singleplotER(cfg,avData(1,1))

FT_TIMELOCKGRANDAVERAGE
FT_MULTIPLOTER
FT_SINGLEPLOTER
