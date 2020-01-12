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


% maybe look for the effect at lateral electrodes for static and central for dynamic? 

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
        end
        subplot(2,3,ff)
        pdPhasePlot(phaseDataToPlot,confRadius);
        %%%% compute difference between left and right adaptation
        % why computing amplitude difference??
        ampDiff(fq,ss,ff) = Axx(conditions(ff,2)).amp(iElec,iFr(fq)) - Axx(conditions(ff,3)).amp(iElec,iFr(fq));
        % phase and amplitude differences 
        phaseDiff(fq,ss,ff) = rad2deg(angle(phaseDataToPlot(2)) - angle(phaseDataToPlot(3)));
        % sin and cos differences to get phiDiff should be equal to
        % phaseDiff... heu pas vraiment... 
        cosDiff = Axx(conditions(ff,2)).cos(iElec,iFr(fq)) - Axx(conditions(ff,3)).cos(iElec,iFr(fq));
        sinDiff = Axx(conditions(ff,2)).sin(iElec,iFr(fq)) - Axx(conditions(ff,3)).sin(iElec,iFr(fq));
        phiDiff(fq,ss,ff) = rad2deg(atan2(cosDiff,sinDiff)); 
    end
%     saveas(gcf,['figures' filesep 'Result_S' num2str(ss) '_Fq' num2str(fq) '.png'],'png')
    end
    dataSbj(:,ss) = Axx;
end
    
% amp(iChan,:) = abs(complexNb);
% cos(iChan,:) = real(complexNb);
% sin(iChan,:) = -imag(complexNb);
    

%%%%%%%%%%%%%%%%%%%%%
% substracting the no-adaptation condition adds noise in the signal so not
% a great thing to do

alpha = deg2rad(squeeze(phaseDiff(1,:,:)));

figure;
for iCond = 1:6
    subplot(2,3,iCond); 
    [phi(iCond) uCI(iCond) lCI(iCond)] = circ_mean(alpha(:,iCond)); % mean direction + upper and lower 95% CI
    circ_plot(alpha(:,iCond),'pretty','ro',true,'linewidth',2,'color','r');
end
saveas(gcf,['figures' filesep 'phaseDiff-adaptation.png'],'png') 
% phi = [2.9314   -2.4366    3.0632   -1.9449   -3.1292   -2.5163]; % mean direction
% upCl = [NaN,-3.34767959928740,NaN,-2.65245642630042,-3.33621068191455,-3.12983270982153];
% confDiff = phi - upCl; % CI for direction 

%%%% comparison phase difference between conditions:
% G1 vs G2h / G1 vs G2l / G2l vs G2h at 10 and 90
% for no adaptation, left, right 
% + between 10 and 90 for G1, G2h G2l

% 10 vs 90 deg test (ie static vs dynamic motion)
% I am looking for a flip in the phase which would correspond to the flip
% in percept
clear dynPhaseDiff
iElec = [22:24; 9:11; 38:40]; iFr = dataSbj(1,1).i1f1;
iCond = [1 10; 4 13; 7 16; 2 3; 5 6; 8 9; 11 12; 14 15; 17 18];
for ee = 1:length(iElec)
for ss=1:size(dataSbj,2)
    for cc=1:length(iCond)
        rad1 = angle(complex(mean(dataSbj(iCond(cc,1),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,1),ss).sin(iElec(ee,:),iFr))));
        rad2 = angle(complex(mean(dataSbj(iCond(cc,2),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,2),ss).sin(iElec(ee,:),iFr))));
        dynPhaseDiff(ss,cc,ee) = rad1-rad2;
    end
end
end
titres = {'noG1','leftG1','rightG1','noG2h','leftG2h','rightG2h','noG2l','leftG2l','rightG2l'};
for ee = 1:length(iElec)
figure;
for iCond = 1:9
    subplot(3,3,iCond); hold on;
    circ_plot(dynPhaseDiff(:,iCond,ee),'pretty','ro',true,'linewidth',2,'color','r');
    title(titres{iCond});
end
saveas(gcf,['figures' filesep 'phaseDiff10vs90-E' num2str(ee) '.png'],'png')
end

% G2l vs G2h 
% here too, some should be the same some should flip if follows the percept
for ee = 1:length(iElec)
for ss=1:size(dataSbj,2)
    iCond = 0;
    for cc=[2 3 5 6 8 9]
        iCond=iCond+1;
        rad1 = angle(complex(mean(dataSbj(cc,ss).cos(iElec(ee,:),iFr)),mean(dataSbj(cc,ss).sin(iElec(ee,:),iFr))));
        rad2 = angle(complex(mean(dataSbj(cc+9,ss).cos(iElec(ee,:),iFr)),mean(dataSbj(cc+9,ss).sin(iElec(ee,:),iFr))));
        phaseLvsH(ss,iCond,ee) = rad1-rad2;
    end
end
end
titres = {'10no','90n','10l','90l','10r','90r'};
for ee = 1:length(iElec)
figure;
for iCond = 1:6
    subplot(3,2,iCond); hold on;
    circ_plot(phaseLvsH(:,iCond,ee),'pretty','ro',true,'linewidth',2,'color','r');
    title(titres{iCond});
end
saveas(gcf,['figures' filesep 'phaseDiffLvsH-E' num2str(ee) '.png'],'png')
end

% Gl vs G2h 
% here too, some should be the same some should flip if follows the percept
iCond = [1 2; 4 5; 7 8; 10 3; 13 6; 16 9];
for ee = 1:length(iElec)
for ss=1:size(dataSbj,2)
    for cc=1:length(iCond)
        rad1 = angle(complex(mean(dataSbj(iCond(cc,1),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,1),ss).sin(iElec(ee,:),iFr))));
        rad2 = angle(complex(mean(dataSbj(iCond(cc,2),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,2),ss).sin(iElec(ee,:),iFr))));
        phase1vsH(ss,cc,ee) = rad1-rad2;
    end
end
end
titres = {'10no','10l','10r','90no','90l','90r',};
for ee = 1:length(iElec)
figure;
for iCond = 1:6
    subplot(2,3,iCond); hold on;
    circ_plot(phase1vsH(:,iCond,ee),'pretty','ro',true,'linewidth',2,'color','r');
    title(titres{iCond});
end
saveas(gcf,['figures' filesep 'phaseDiff1vsH-E' num2str(ee) '.png'],'png')
end

% Gl vs G2l 
% here too, some should be the same some should flip if follows the percept
iCond = [1 11; 4 14; 7 17; 10 12; 13 15; 16 18];
for ee = 1:length(iElec)
for ss=1:size(dataSbj,2)
    for cc=1:length(iCond)
        rad1 = angle(complex(mean(dataSbj(iCond(cc,1),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,1),ss).sin(iElec(ee,:),iFr))));
        rad2 = angle(complex(mean(dataSbj(iCond(cc,2),ss).cos(iElec(ee,:),iFr)),mean(dataSbj(iCond(cc,2),ss).sin(iElec(ee,:),iFr))));
        phase1vsL(ss,cc,ee) = rad1-rad2;
    end
end
end
titres = {'10no','10l','10r','90no','90l','90r',};
for ee = 1:length(iElec)
figure;
for iCond = 1:6
    subplot(2,3,iCond); hold on;
    circ_plot(phase1vsL(:,iCond,ee),'pretty','ro',true,'linewidth',2,'color','r');
    title(titres{iCond});
end
saveas(gcf,['figures' filesep 'phaseDiff1vsL-E' num2str(ee) '.png'],'png')
end