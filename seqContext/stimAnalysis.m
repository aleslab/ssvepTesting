
clearvars;
dataIn = '/Users/marleneponcet/Documents/data/seqContext/cleanData/';
filesDir = 'C:\Users\Marlene\Documents\dataStAndrews\seqContext\clean\';
filesList = dir([filesDir '*stim.mat']);

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

figure; hold on
for condition = 1:length(allCond)
    plot(avData(ss,condition).time,avData(ss,condition).avg(23,:),'LineWidth',2)
end
legend('motion','seq','reversal','two','rand')

% figure; hold on
% for condition = 1:length(allCond)
%     plot(avDataT(23,:,condition),'LineWidth',2)
% end

end

for condition = 1:length(allCond)
    grandavg(condition) = ft_timelockgrandaverage(cfg,avData(1,condition),avData(2,condition),avData(3,condition));
end
figure; hold on
for condition = 1:length(allCond)
    plot(grandavg(condition).time,grandavg(condition).avg(23,:),'LineWidth',2)
end
legend('motion','seq','reversal','two','rand')

ft_singleplotER(cfg,avData(1,1))

FT_TIMELOCKGRANDAVERAGE
FT_MULTIPLOTER
FT_SINGLEPLOTER
