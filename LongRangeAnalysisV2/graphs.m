condition = 10;
channel = 23;

% time-frequency plot
plot(squeeze(gpData(condition).freq),squeeze(gpData(condition).amp(channel,:)))
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)))

figure;
condition = 9;
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)),0.8)
hold on;
condition = 10;
bar(squeeze(gpData(condition).freq((gpData(condition).freq<30 & gpData(condition).freq>0))),squeeze(gpData(condition).amp(channel,2:60)),0.3,'FaceColor','r')

ft_topoplotTFR(cfg,squeeze(gpData(condition).amp(2:end,10)))
plotTopo(squeeze(gpData(condition).amp(:,10)),cfg.layout)