addpath C:\Users\Marlene\Documents\git\fieldtrip-aleslab-fork
addpath C:\Users\Marlene\Documents\git\ssvepTesting\svndlCopy
addpath C:\Users\Marlene\Documents\git\ssvepTesting\biosemiUpdated
addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFunctions
addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFunctions\CircStat2012a
ft_defaults
dataDir = 'C:\Users\Marlene\Documents\dataStAndrews\MAE\AxxPooled\';


addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions/CircStat2012a
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
    iFr   = [Axx(1).i1f1 Axx(1).i1f1*2-1];
    conditions = [1 4 7; 2 5 8; 11 14 17; 10 13 16; 3 6 9; 12 15 18];
    for fq=1:length(iFr)
        figure;hold on;
    for ff = 1:6 % for each test stimulus
        clear confRadius phaseDataToPlot
        for iCond = 1:size(conditions,2) % 3 adaptations
            phaseDataToPlot(iCond) = complex( Axx(conditions(ff,iCond)).cos(iElec,iFr(fq)), Axx(conditions(ff,iCond)).sin(iElec,iFr(fq)));
            %TODO: Fix this! this is an underestimate of the 95% confidence intervals!
            confRadius(iCond) = Axx(conditions(ff,iCond)).confradius(iElec,iFr(fq));
            alpha(ss,fq,ff,iCond) = angle(phaseDataToPlot(iCond));
       end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
    end
%     saveas(gcf,['figures' filesep 'Result_S' num2str(ss) '_Fq' num2str(fq) '.png'],'png')
    end
    dataSbj(:,ss) = Axx;
end

%%%%%%%%%%%%%%%%%%
%%% average
[groupAv,sbjProj] = averageAxxWithStd(dataSbj);
interactiveSteadyStatePlot2(cfg,groupAv)

%%% figure with all conditions (no-left-right adaptation) on the same
%%% phasor
iFr   = [groupAv(1).i1f1 groupAv(1).i1f1*2-1];
conditions = [1 4 7; 2 5 8; 11 14 17; 10 13 16; 3 6 9; 12 15 18];
listTitles = {'G1-10','G2-10h','G2-10l','G1-90','G2-90h','G2-90l'};
for fq=1:length(iFr)
for iElec = [10 23 39]
    figure;hold on;
    for ff = 1:6
        clear confRadius phaseDataToPlot
        for iCond = 1:size(conditions,2)
            phaseDataToPlot(iCond) = complex( groupAv(conditions(ff,iCond)).cos(iElec,iFr(fq)), groupAv(conditions(ff,iCond)).sin(iElec,iFr(fq)));
            %TODO: Fix this! this is an underestimate of the 95% confidence intervals!
            confRadius(iCond) = groupAv(conditions(ff,iCond)).confradius(iElec,iFr(fq));
        end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
        title(listTitles{ff})
    end
    % colours: red=left, yellow=right, blue=none
    saveas(gcf,['figures' filesep 'Result_Average_Fq' num2str(fq) '_E' num2str(iElec) '.png'],'png')
end
end

%%% topographies
% C1 = 1 grating, C2 = 2 gratings high SF, C3 = 2 gratings low SF
% conditions: none left right for top10 bottom90 
nameCond = {'0.5_SF', '0.5+1_SF','0.5+0.25_SF'};
conditions = [1 4 7 10 13 16; 2 5 8 3 6 9; 11 14 17 12 15 18];
for fq=1:length(iFr)
for ff=1:size(conditions,1)
    figure;hold on;
    for iCond = 1:length(conditions)
        subplot(2,3,iCond)
        plotTopo(groupAv(conditions(ff,iCond)).amp(:,iFr(fq)),cfg.layout);
        colorbar
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 1500 1000])
    colormap('hot');
    saveas(gcf,['figures' filesep 'Result_Topo_' nameCond{ff} 'fq' num2str(fq) '.png'],'png')
end
end
    
%%%
%%%% gonna just plot the ERP for now... 
condName = {'G1-10','G1-90','G2h-10','G2h-90','G2l-10','G2l-90'};
iElec = [22:24; 9:11; 38:40];
for ee=1:length(iElec)
    ii=1;
    for iCond= [1 10 2 3 11 12]
        figure; hold on;
        plot(groupAv(iCond).time,mean(groupAv(iCond).wave(iElec(ee,:),:),1),'LineWidth',1)
        plot(groupAv(iCond+3).time,mean(groupAv(iCond+3).wave(iElec(ee,:),:),1),'LineWidth',1)
        plot(groupAv(iCond+6).time,mean(groupAv(iCond+6).wave(iElec(ee,:),:),1),'LineWidth',1)
        legend('none','left','right');
        title(condName{ii});
        saveas(gcf,['figures' filesep 'Result_ERP_' condName{ii} 'elec' num2str(ee) '.png'],'png');
        ii = ii+1;
    end
end





%%%%%% try to get variance for the phase amplitude using sbjProj and 
%%%%%% the variance for the direction using circStat 


elec = 23; freq = groupAv(1).i1f1;
% sbjProj(elec,freq,iCond,ss)

% get the mean CI direction
% use variable alpha
% each figure plot each condition, no adapt, left, right
for fq=1:2
for ff=1:6
     figure; 
for iCond = 1:3
    subplot(1,3,iCond); hold on;
    [phi(fq,ff,iCond) uCI(fq,ff,iCond) lCI(fq,ff,iCond)] = circ_mean(alpha(:,fq,ff,iCond)); % mean direction + upper and lower 95% CI
    circ_plot(alpha(:,fq,ff,iCond),'pretty','ro',true,'linewidth',2,'color','r');
end
end
end

cc=1;
for ff = 1:6
    for iCond=1:3
        stdCond(cc) = std(sbjProj(elec,freq,cc,:));
        cc=cc+1;
    end
end
% mean direction: phi
% 95CI: uCI
% mean amplitude: groupAv.amp
% stdCond
% test C1, no adapt
compMean = complex(groupAv(1).amp(elec,freq), phi(1,1,1));
%%%%%% NO! amplitude is not cosine!!!
confMean = complex(stdCond(1), uCI(1,1,1));
figure;pdPhasePlot(compMean,confMean);

cosMean = cos(phi(1,1,1)) * groupAv(1).amp(elec,freq);
compMean = complex(0.9269, phi(1,1,1));


% convert angles to unit vectors
z = exp(1i*phi(1,1,1));
figure;pdPhasePlot(z,confMean);



% amp(iChan,:) = abs(complexNb);
% cos(iChan,:) = real(complexNb);
% sin(iChan,:) = -imag(complexNb);




% %%%% from mrCurrent code
% tNellip = 30;
% tTh = linspace( 0, 2*pi, tNellip )';
% tNSbjs = length(listData);
% if SEM
%     tNormK = 1 / (tNSbjs-2);
% else % 95% CI
%     tNormK = (tNSbjs-1)/tNSbjs/(tNSbjs-2) * finv( 0.95, 2, tNSbjs - 2 );
% end
% [ tEVec, tEVal ] = eig( cov( [ real( tYSubj ), imag( tYSubj ) ] ) );	% Compute eigen-stuff
% tXY = [ cos(tTh), sin(tTh) ] * sqrt( tNormK * tEVal ) * tEVec';		% Error/confidence ellipse
                                

% average across all odd/even harmonics
%%%% BUT THIS IS NOT POSSIBLE FOR THE PHASE!!! Phase might be different at
%%%% different harmonics so cannot just average them! Can average the
%%%% difference though






% I wonder, could I just fake one condition to be at 0deg phase so that all
% subjects are aligned on the same phase? e.g. left adaptation is at 45deg
% for S01, for all conditions I remove 45deg. then expect right adpatation
% to be at 180deg.

% how to pick electrodes / frequencies?
% use average of a few electrodes instead of just 23? 
% average of 3 electrodes does not change much
%%% phasor
iFr   = [groupAv(1).i1f1 groupAv(1).i1f1*2-1];
conditions = [1 4 7; 2 5 8; 11 14 17; 10 13 16; 3 6 9; 12 15 18];
listTitles = {'G1-10','G2-10h','G2-10l','G1-90','G2-90h','G2-90l'};
iElec=[22:24; 9:11; 38:40];
for fq=1:length(iFr)
for iE = 1:3
    figure;hold on;
    for ff = 1:6
        clear confRadius phaseDataToPlot
        for iCond = 1:size(conditions,2)
            phaseDataToPlot(iCond) = complex( mean(groupAv(conditions(ff,iCond)).cos(iElec(iE,:),iFr(fq))), mean(groupAv(conditions(ff,iCond)).sin(iElec(iE,:),iFr(fq))));
            %TODO: Fix this! this is an underestimate of the 95% confidence intervals!
            confRadius(iCond) = NaN;
        end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
        title(listTitles{ff})
    end
    % colours: red=left, yellow=right, blue=none
    saveas(gcf,['figures' filesep 'Result_Average_Fq' num2str(fq) '_meanE' num2str(mean(iElec)) '.png'],'png')
end
end






