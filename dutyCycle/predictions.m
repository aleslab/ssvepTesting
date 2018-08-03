clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = 1:length(listData); 
keepSbj = [1:5 7 9:11 13:15]; 


% compute average (for later comparison)
for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    dataSbj(:,ss) = Axx;
end
avData = averageAxx(dataSbj(:,keepSbj));    
% % compute some stats
% test = gpStatT2_DC(dataSbj,avData);

% compute predictions (where possible)
for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);

    % use cond 11 = 2Hz 12.5% to predict cond 7 = 5Hz 25% 
    tmp = sumAxxWithShift(Axx(11),Axx(11),length(Axx(7).time)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred ' Axx(7).condLabel];
    pred(1,ss) = tmp; % need temporary variable until I get the same nb of fields
    difftmp = computeDiff(tmp,Axx(7)); % compute the difference (reconstructed always first)
    difftmp.condLabel = ['diff ' Axx(7).condLabel];
    diffmotion(1,ss) = difftmp;
    clear tmp;difftmp;
    
    % use cond 12 = 2Hz 25% to predict cond 8 = 5Hz 50% 
    tmp = sumAxxWithShift(Axx(12),Axx(12),length(Axx(8).time)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred ' Axx(8).condLabel];
    pred(2,ss) = tmp; % need temporary variable until I get the same nb of fields
    difftmp = computeDiff(tmp,Axx(8)); % compute the difference (reconstructed always first)
    difftmp.condLabel = ['diff ' Axx(8).condLabel];
    diffmotion(2,ss) = difftmp;
    clear tmp;difftmp;
    
    % use cond 7 = 5Hz 25% to predict cond 3 = 10Hz 50% 
    tmp = sumAxxWithShift(Axx(7),Axx(7),length(Axx(3).time)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred ' Axx(3).condLabel];
    pred(3,ss) = tmp; % need temporary variable until I get the same nb of fields
    difftmp = computeDiff(tmp,Axx(3)); % compute the difference (reconstructed always first)
    difftmp.condLabel = ['diff ' Axx(3).condLabel];
    diffmotion(3,ss) = difftmp;
    clear tmp;difftmp;

    % use cond 6 = 5Hz 12.5% to predict cond 2 = 10Hz 25% 
    tmp = sumAxxWithShift(Axx(6),Axx(6),length(Axx(2).time)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred ' Axx(2).condLabel];
    pred(4,ss) = tmp; % need temporary variable until I get the same nb of fields
    difftmp = computeDiff(tmp,Axx(2)); % compute the difference (reconstructed always first)
    difftmp.condLabel = ['diff ' Axx(2).condLabel];
    diffmotion(4,ss) = difftmp;
    clear tmp;difftmp;
    
    % use cond 11 = 2Hz 12.5% to predict cond 3 = 10Hz 50% 
    tmp = sumAxxWithShift(Axx(11),sumAxxWithShift(Axx(11),Axx(11)),length(Axx(3).time)); 
    tmp.condLabel = ['pred2 ' Axx(3).condLabel];
    pred(5,ss) = tmp;    
    clear tmp;

    %%%%%%%%% facking motion (not spatial difference)
     % use cond 11 = 2Hz 12.5% to predict cond 17  
    tmp = sumAxxWithShift(Axx(11),Axx(11)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred' Axx(17).condLabel];
    pred(6,ss) = tmp; % need temporary variable until I get the same nb of fields
%     difftmp = computeDiff(tmp,Axx(17)); % compute the difference (reconstructed always first)
%     difftmp.condLabel = ['diff ' Axx(17).condLabel];
%     diffmotion(1,ss) = difftmp;
    clear tmp difftmp;
    
    % use cond 12 = 2Hz 25% to predict cond 22 
    tmp = sumAxxWithShift(Axx(12),Axx(12)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred' Axx(22).condLabel];
    pred(7,ss) = tmp; % need temporary variable until I get the same nb of fields
%     difftmp = computeDiff(tmp,Axx(22)); % compute the difference (reconstructed always first)
%     difftmp.condLabel = ['diff ' Axx(22).condLabel];
%     diffmotion(2,ss) = difftmp;
    clear tmp difftmp;
    
%     % use cond 13 = 2Hz 50% to predict cond 19 
%     tmp = sumAxxWithShift(Axx(13),Axx(13)); % 3rd argument optional = size of the new cycle
%     tmp.condLabel = ['pred' Axx(19).condLabel];
%     pred(8,ss) = tmp; % need temporary variable until I get the same nb of fields
% %     difftmp = computeDiff(tmp,Axx(19)); % compute the difference (reconstructed always first)
% %     difftmp.condLabel = ['diff ' Axx(19).condLabel];
% %     diffmotion(3,ss) = difftmp;
%     clear tmp difftmp;
    
    % use cond 7 = 5Hz 25% to predict cond 18
    tmp = sumAxxWithShift(Axx(7),Axx(7)); % 3rd argument optional = size of the new cycle
    tmp.condLabel = ['pred' Axx(18).condLabel];
    pred(8,ss) = tmp; % need temporary variable until I get the same nb of fields
%     difftmp = computeDiff(tmp,Axx(18)); % compute the difference (reconstructed always first)
%     difftmp.condLabel = ['diff ' Axx(18).condLabel];
%     diffmotion(4,ss) = difftmp;
    clear tmp difftmp;
end

avPred = averageAxx(pred(:,keepSbj));
avDiffMotion = averageAxx(diffmotion(:,keepSbj));
interactiveSteadyStatePlot2(cfg,[avData,avPred,avDiffMotion])






%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [7 8 3 2];channel=23;
position=[3 4 2 1];

figure;
iff=2;
figure;
for cond=1:4
    subplot(2,2,position(cond)); hold on;
    plotTopo(avDiffMotion(cond).amp(:,(avDiffMotion(cond).i1f1-1)*iff+1),cfg.layout);
    colorbar;
    caxis([0 1.6]);
end
saveas(gcf,'predDiffTopo.jpg')

channel = [9 39];
for chan=1:2
figure;
for cond=1:4
    subplot(2,2,position(cond)); hold on;
    bar(squeeze(avDiffMotion(cond).freq((avDiffMotion(cond).freq<15.5 & avDiffMotion(1).freq>0))),squeeze(avDiffMotion(cond).amp(channel(chan),2:36)))
    ylim([0 2]);
end
saveas(gcf,['predDiffFreq' num2str(channel(chan)) '.jpg'])
end