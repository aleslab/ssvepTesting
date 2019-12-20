
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/MAEv3/pooled/';

listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};



for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    %%%% results per participant
    iElec = 23;
    iFr   = Axx(1).i1f1;
    figure;hold on;
    conditions = [1 4 7; 2 5 8; 11 14 17; 10 13 16; 3 6 9; 12 15 18];
    for ff = 1:6
        clear confRadius phaseDataToPlot
        for iCond = 1:size(conditions,2)
            phaseDataToPlot(iCond) = complex( Axx(conditions(ff,iCond)).cos(iElec,iFr), Axx(conditions(ff,iCond)).sin(iElec,iFr));
            %TODO: Fix this! this is an underestimate of the 95% confidence intervals!
            confRadius(iCond) = Axx(conditions(ff,iCond)).confradius(iElec,iFr);
        end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
    end
    saveas(gcf,['figures' filesep 'Result_S' num2str(ss) '.png'],'png')
    
    dataSbj(:,ss) = Axx;
end

%%%%%%%%%%%%%%%%%%
%%% average
[groupAv,sbjProj] = averageAxxWithStd(dataSbj);
interactiveSteadyStatePlot2(cfg,groupAv)

%%% figure with all conditions (no-left-right adaptation) on the same
%%% phasor
iFr   = groupAv(1).i1f1;
conditions = [1 4 7; 2 5 8; 11 14 17; 10 13 16; 3 6 9; 12 15 18];
listTitles = {'G1-10','G2-10h','G2-10l','G1-90','G2-90h','G2-90l'};
for iElec = [10 23 39]
    figure;hold on;
    for ff = 1:6
        clear confRadius phaseDataToPlot
        for iCond = 1:size(conditions,2)
            phaseDataToPlot(iCond) = complex( groupAv(conditions(ff,iCond)).cos(iElec,iFr), groupAv(conditions(ff,iCond)).sin(iElec,iFr));
            %TODO: Fix this! this is an underestimate of the 95% confidence intervals!
            confRadius(iCond) = groupAv(conditions(ff,iCond)).confradius(iElec,iFr);
        end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
        title(listTitles{ff})
    end
    % colours: red=left, yellow=right, blue=none
    saveas(gcf,['figures' filesep 'Result_Average_E' num2str(iElec) '.png'],'png')
end

%%% topographies
conditions = [1 4 7 10 13 16; 2 5 8 3 6 9; 11 14 17 12 15 18];
for ff=1:size(conditions,1)
    figure;hold on;
    for iCond = 1:length(conditions)
        subplot(2,3,iCond)
        plotTopo(groupAv(conditions(ff,iCond)).amp(:,iFr),cfg.layout);
        colorbar
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1500 1000])
    colormap('hot');
    saveas(gcf,['figures' filesep 'Result_Topo_Cond' num2str(ff) '.png'],'png')
end
% conditions: none left right for top10 bottom90 
% C1=G1, C2=G2 high, C3=G2low


    
% average across all odd/even harmonics


%%%%%%%%%%%%%%%%%%%%%
%%% substraction left-right
% 2 ways should give the same results

% phase and amplitude difference 
% average

% sin and cos difference
% average

%%%%%%%%%%%%%%%%%%%%%
% substracting the no-adaptation condition adds noise in the signal so not
% a great thing to do




% I wonder, could I just fake one condition to be at 0deg phase so that all
% subjects are aligned on the same phase? e.g. left adaptation is at 45deg
% for S01, for all conditions I remove 45deg. then expect right adpatation
% to be at 180deg.

% how to pick electrodes / frequencies?

