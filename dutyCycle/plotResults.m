
clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};


for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    dataSbj(:,ss) = Axx;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SNR for each participant (check if needs to be rejected)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum all the harmonics and divide by the neighbouring harmonics
% doesn't give the same results if do the snr per cond then average across
% cond
clear snr;
for ss=1:size(dataSbj,2)
    signal=0;noise=0;
    for cond = 1:15
        nbHarm = 0;
        currentHarm=1;
        currentFq = (dataSbj(cond,ss).i1f1-1)*currentHarm+1;
        while dataSbj(cond,ss).freq(:,currentFq) < 85
            signal = signal + dataSbj(cond,ss).amp(:,currentFq);
            noise = noise + ((dataSbj(cond,ss).amp(:,currentFq-1) + dataSbj(cond,ss).amp(:,currentFq+1)) /2);
            currentHarm = currentHarm+1;
            currentFq = (dataSbj(cond,ss).i1f1-1)*currentHarm+1;
        end
    end
    snr(:,ss) = signal./noise;
end

figure
for ss=1:size(dataSbj,2)
subplot(3,5,ss); hold on;
plotTopo(snr(:,ss),cfg.layout);
caxis([0.2 2.2]);
end
colorbar
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'snrSbj.jpg')


%% average
keepSbj = [1:size(dataSbj,2)]; % reject S6?
keepSbj = [1:5 7 9:11 13]; 

avData = averageAxx(dataSbj(:,keepSbj));

% compute some stats
% avData = gpStatT2_DC(dataSbj,avData);

interactiveSteadyStatePlot2(cfg,avData) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot topographies for different conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% f1
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,avData(cond).i1f1),cfg.layout);
    caxis([0 2.5]);
end
colorbar; 
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'f1topo.jpg')
%% f2
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*2+1),cfg.layout);
    caxis([0 1.8]); 
end
colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'f2topo.jpg')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%% f3
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*3+1),cfg.layout);
    caxis([0 1]); 
end
colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'f3topo.jpg')
%% f4
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*4+1),cfg.layout);
    caxis([0 1]); 
end
colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'f4topo.jpg')
%% f5
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*5+1),cfg.layout);
    caxis([0 0.5]); 
end
colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'f5topo.jpg')


%%%%%%%%%%%%%%%%%% Plot topographies for MOTION
position = [11 7 3 9 15 8 13];
for iff=1:5
    figure;
    for cond=16:22
        subplot(3,5,position(cond-15)); hold on;
        plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*iff+1),cfg.layout);
        %     caxis([0 2.5]);
        colorbar;
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 2000 2500])
    saveas(gcf,[num2str(iff) 'f1topoMotion.jpg'])
end



%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Electrode based plots
%%%%%%%%%%%%%%%%%%

% plot mean amplitude for the 5 duty cycle at the 5 harmonics
% in 3 graphs separately for the 3 frequencies

pickElec = 23; % 23 Oz, 9, B6=38
for cond=1:15 % exclude motion conditions
    for ff = 1 : 5
%         avData(cond).freq(avData(cond).i1f1*ff)
        meanAmp(cond,ff)=avData(cond).amp(pickElec,(avData(cond).i1f1-1)*ff+1);
    end
end
% figure;bar(meanAmp(:,:));

%%
figure;
for ff=1:3
    subplot(1,3,ff)
    bar(transpose([meanAmp(1:5,ff),meanAmp(6:10,ff),meanAmp(11:15,ff)]))
    ylim([0 2.5])
    xticklabels({'10.6Hz','5.3Hz','2.7Hz'})
    legend('12.5% DC','25% DC','50% DC','75% DC','87.5% DC','Location','NorthWest')
    title(['mean amplitude ' num2str(ff) 'f1'])
    xlabel('presentation rate')
end
saveas(gcf,[num2str(pickElec) ' Amplitude.jpg'])

figure;
for ff=1:4
subplot(2,2,ff)
imagesc(transpose([meanAmp(1:5,ff),meanAmp(6:10,ff),meanAmp(11:15,ff)]));
colorbar;caxis([0 2.5]); 
xticks(1:5)
xticklabels({'12.5% DC','25% DC','50% DC','75% DC','87.5% DC'})
yticks(1:3)
yticklabels({'10.6Hz','5.3Hz','2.7Hz'})
title([num2str(pickElec) ' amplitude ' num2str(ff) 'f1'])
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [200 200 1000 800])
saveas(gcf,[num2str(pickElec) 'AmplitudeConfMat.jpg'])

%%
%%%%%%%%%%%%%%%%%% plot amplitude at one electrode for different conditions
%%%%%%%%%%%%%%%%%% depending on timeON, timeOff, duty cycle
vect = [1 2 4 6 7];
figure;line(1/85*vect*1000,meanAmp(1:5,1),'Color','b','Marker','o')
hold on; line(1/85*2*vect*1000,meanAmp(6:10,1),'Color','r','Marker','o')
hold on; line(1/85*4*vect*1000,meanAmp(11:15,1),'Color','g','Marker','o')
line(1/85*vect*1000,meanAmp(1:5,2),'LineStyle','--','Color','b','Marker','*')
hold on; line(1/85*2*vect*1000,meanAmp(6:10,2),'LineStyle','--','Color','r','Marker','*')
hold on; line(1/85*4*vect*1000,meanAmp(11:15,2),'LineStyle','--','Color','g','Marker','*')
xlabel('time ON')
ylabel('amplitude')
legend('10Hz 1f1','5Hz 1f1','2Hz 1f1','10Hz 2f1','5Hz 2f1','2Hz 2f1')
title('Time On')
saveas(gcf,[num2str(pickElec) ' TimeON.jpg'])


figure;line(1/85.*(fliplr(vect))*1000,meanAmp(1:5,1),'Color','b','Marker','o')
hold on; line(1/85*(fliplr(2*vect))*1000,meanAmp(6:10,1),'Color','r','Marker','o')
line(1/85*(fliplr(4*vect))*1000,meanAmp(11:15,1),'Color','g','Marker','o')
line(1/85*fliplr(vect)*1000,meanAmp(1:5,2),'LineStyle','--','Color','b','Marker','*')
line(1/85*fliplr(2*vect)*1000,meanAmp(6:10,2),'LineStyle','--','Color','r','Marker','*')
line(1/85*fliplr(4*vect)*1000,meanAmp(11:15,2),'LineStyle','--','Color','g','Marker','*')
xlabel('time OFF')
ylabel('amplitude')
legend('10Hz 1f1','5Hz 1f1','2Hz 1f1','10Hz 2f1','5Hz 2f1','2Hz 2f1')
title('Time Off')
saveas(gcf,[num2str(pickElec) 'TimeOff.jpg'])

figure;line((1:5)/8*100,meanAmp(1:5,1),'Color','b','Marker','o')
hold on; line((1:5)/8*100,meanAmp(6:10,1),'Color','r','Marker','o')
line((1:5)/8*100,meanAmp(11:15,1),'Color','g','Marker','o')
line((1:5)/8*100,meanAmp(1:5,2),'LineStyle','--','Color','b','Marker','*')
line((1:5)/8*100,meanAmp(6:10,2),'LineStyle','--','Color','r','Marker','*')
line((1:5)/8*100,meanAmp(11:15,2),'LineStyle','--','Color','g','Marker','*')
legend('10Hz 1f1','5Hz 1f1','2Hz 1f1','10Hz 2f1','5Hz 2f1','2Hz 2f1')
xlabel('Duty Cycle ratio')
ylabel('amplitude')
title('Duty cycle')
saveas(gcf,[num2str(pickElec) 'DutyCycle.jpg'])




%% Mean total signal per condition
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% mean total signal (density) = sqrt(total signal power) / nb of harmonics
% total sig power = f1*f1 + 2f1*2f1 + 3f1*3f1 + ....
% signal LP at 85 Hz, so consider freq < 85

for cond = 1:22
sigPower = zeros;
nbHarm = 0;
currentHarm=1;
currentFq = (avData(cond).i1f1-1)*currentHarm+1;
while avData(cond).freq(:,currentFq) < 85
    nbHarm = nbHarm +1;
%     avData(cond).freq(:,currentFq)
    sigPower = sigPower + avData(cond).amp(:,currentFq).*avData(cond).amp(:,currentFq);
    currentHarm = currentHarm+1; 
    currentFq = (avData(cond).i1f1-1)*currentHarm+1;
end
power(:,cond) = sigPower;
meanTotSig(:,cond) = sqrt(sigPower) / nbHarm;
end

figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(power(:,cond),cfg.layout);
    colorbar; 
%     caxis([0 0.3]);
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'sigPowerTopo.jpg')

figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(meanTotSig(:,cond),cfg.layout);
%     caxis([0 0.3]);
colorbar;
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'meanTotSigTopo.jpg')



%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% plot MOTION
%%%%%%%%%%%%%%%%%%%%
position = [11 7 3 9 15 8 13];
for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plotTopo(meanTotSig(:,cond+15),cfg.layout);
%     caxis([0 0.3]);
colorbar;
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'MotionMeanTotSigTopo.jpg')

for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plotTopo(power(:,cond+15),cfg.layout);
%     caxis([0 0.3]);
colorbar;
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'MotionSigPowerTopo.jpg')














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% waveforms and phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pickElec = 23; % 23 Oz, 9, B6=38
figure;hold on;
subplot(3,1,1);hold on;
for cond = 1:5
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'LineWidth',2)
end
legend('12','25','50','75','87')
title('10Hz')
subplot(3,1,2);hold on;
for cond = 6:10
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'LineWidth',2)
end
legend('12','25','50','75','87')
title('5Hz')
subplot(3,1,3);hold on;
for cond = 11:15
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'LineWidth',2)
end
legend('12','25','50','75','87')
title('2Hz')
saveas(gcf,'waveformOz.jpg')

