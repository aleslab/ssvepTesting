clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7 9:11 13:15]; 


% compute average (for later comparison)
for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    dataSbj(:,ss) = Axx;
end
avData = averageAxx(dataSbj(:,keepSbj));    

clear diff

tmpAxx = avData(16);
tmpAxx.wave = avData(16).wave(:,1:size(avData(16).wave,2)/2);
diff(1) = computeDiff(avData(11),tmpAxx);
% diff(1).condLabel = ['diff' Axx(11).condLabel];

tmpAxx = avData(17);
tmpAxx.wave = avData(17).wave(:,1:size(avData(17).wave,2)/2);
diff(2) = computeDiff(avData(7),tmpAxx);
% diff(2).condLabel = ['diff' Axx(7).condLabel];

tmpAxx = avData(18);
tmpAxx.wave = avData(18).wave(:,1:size(avData(18).wave,2)/2);
diff(3) = computeDiff(avData(3),tmpAxx);
% diff(3).condLabel = ['diff' Axx(3).condLabel];

tmpAxx = avData(19);
tmpAxx.wave = avData(19).wave(:,1:size(avData(19).wave,2)/2);
diff(4) = computeDiff(avData(9),tmpAxx);
% diff(4).condLabel = ['diff' Axx(9).condLabel];

tmpAxx = avData(20);
tmpAxx.wave = avData(20).wave(:,1:size(avData(20).wave,2)/2);
diff(5) = computeDiff(avData(15),tmpAxx);
% diff(5).condLabel = ['diff' Axx(15).condLabel];

tmpAxx = avData(21);
tmpAxx.wave = avData(21).wave(:,1:size(avData(21).wave,2)/2);
diff(6) = computeDiff(avData(13),tmpAxx);
% diff(6).condLabel = ['diff' Axx(13).condLabel];

tmpAxx = avData(22);
tmpAxx.wave = avData(22).wave(:,1:size(avData(22).wave,2)/2);
diff(7) = computeDiff(avData(8),tmpAxx);
% diff(7).condLabel = ['diff' Axx(8).condLabel];

ttavData = rmfield(avData,'condLabel');
interactiveSteadyStatePlot2(cfg,[ttavData diff])


figure;
position = [11 7 3 9 15 13 8];
for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plotTopo(diff(cond).amp(:,avData(cond).i1f1),cfg.layout);
    caxis([0 0.5]);
% colorbar;
end
colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'diffStatMot.jpg')

figure;
position = [11 7 3 9 15 13 8];
for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plot(diff(cond).wave(23,:));
    ylim([-1 1.5])
    xlim([0 200])
end
saveas(gcf,'diffStatMotwave.jpg')

