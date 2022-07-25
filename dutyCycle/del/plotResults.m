
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
for ss = 1:size(dataSbj,2)
    signal=0;noise=0;
    for cond = 1:length(dataSbj)
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
colorbar; 
% caxis([1 3.5]);
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'snrSbj.jpg')


%% average
keepSbj = [1:size(dataSbj,2)]; % reject S6?
keepSbj = [1:5 7:11 13:15]; 

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
    % caxis([0 2.5]);
    if cond<6
        caxis([0 1.2]);
    else
        caxis([0 2.5]);
    end
    if cond == 5 || cond == 10 || cond == 15
        colorbar;
    end
end
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
position = [11 7 3 9 15 13 8];
limAx = [0.5 2.5 0.5 2 0.5 1];
for iff=1:6
    figure;
    for cond=16:22
        subplot(3,5,position(cond-15)); hold on;
        plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*iff+1),cfg.layout);
            caxis([0 limAx(iff)]);
%         colorbar;
    end
    colorbar
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 2000 2500])
    saveas(gcf,[num2str(iff) 'f1topoMotion.jpg'])
end

% comparison f1 static
    figure;
    for cond=16:22
        subplot(3,5,position(cond-15)); hold on;
        plotTopo(avData(cond).amp(:,(avData(cond).i1f1-1)*2+1),cfg.layout);
    if cond==18
        caxis([0 1.2]);
    else
        caxis([0 2.5]);
    end
    end
    set(gcf,'PaperPositionMode','auto')
    set(gcf, 'Position', [0 0 2000 2500])
    saveas(gcf,'comptopoMotion1f1.jpg')





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

pickElec = 23; % 23 Oz, 9, B6=38, 16=best SNR
figure;hold on;
for cond=1:5:15 % exclude motion conditions
   plot(meanTotSig(pickElec,cond:cond+4));
end
legend('10','5','2');
title('meanTotSig Oz');
figure;hold on;
for cond=1:5:15 % exclude motion conditions
   plot(power(pickElec,cond:cond+4));
end
legend('10','5','2');
title('power Oz');

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

figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(power(:,cond),cfg.layout);
    colorbar; 
%     caxis([0 8]);
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'sigPowerTopo.jpg')

figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    plotTopo(meanTotSig(:,cond),cfg.layout);
    if cond>10
        caxis([0 0.1]);
    else
        caxis([0 0.16]);
    end
    if cond == 5 || cond == 10 || cond == 15
        colorbar;
    end
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'meanTotSigTopo.jpg')



%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% plot MOTION
%%%%%%%%%%%%%%%%%%%%
figure;
position = [11 7 3 9 15 13 8];
for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plotTopo(meanTotSig(:,cond+15),cfg.layout);
    if cond>10
        caxis([0 0.1]);
    else
        caxis([0 0.16]);
    end
    if cond == 5 || cond == 10 || cond == 15
        colorbar;
    end
%     caxis([0 0.3]);
colorbar;
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'MotionMeanTotSigTopo.jpg')

figure;
for cond=1:length(position)
    subplot(3,5,position(cond)); hold on;
    plotTopo(power(:,cond+15),cfg.layout);
%     caxis([0 0.3]);
colorbar;
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'MotionSigPowerTopo.jpg')



%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Correl
%%%%%%%%%%%%%%%%%%
for cond=1:15 %
   avAmp(cond)=avData(cond).amp(23,(avData(cond).i1f1));
end
for cond=16:22 % for motion take the 2f1
   avAmp(cond)=avData(cond).amp(23,(avData(cond).i1f1*2-1));
end
figure;hold on;
scatter(avAmp(1:5),static(1,:),'filled')
scatter(avAmp(6:10),static(2,:),'filled')
scatter(avAmp(11:15),static(3,:),'g','filled')
scatter(avAmp(16:22),moving([3 5 7 11 15 9 8]),'k','filled')
ylabel('motion rating')
xlabel('amplitude')
legend('10','5','2','moving')
saveas(gcf,'ampByrating.jpg')


figure;hold on; 
for flick=1:3
x = static(flick,:);
ss = flick+4*(flick-1);
y = avAmp(ss:ss+4);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'o',x,f,'-') 
end
x = moving([3 5 7 11 15 9 8]);
y = avAmp(16:22);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'o',x,f,'-') 
xlabel('motion rating')
ylabel('amplitude')

figure;hold on; 
y = [static(2,:) static(3,:)];
x = avAmp(6:15);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'o',x,f,'-') 
y = moving([3 5 7 11 15 9 8]);
x = avAmp(16:22);
p = polyfit(x,y,1); 
f = polyval(p,x); 
plot(x,y,'o',x,f,'-') 
ylim([0 3])
legend('static','','moving')
ylabel('motion rating')
xlabel('amplitude')


figure;hold on;
plot(avAmp(1:5),'x--','MarkerSize',10);plot(avAmp(6:10),'x--','MarkerSize',10);plot(avAmp(11:15),'x--','MarkerSize',10)
plot([2 3 4],avAmp([17 22 19]),'.-','MarkerSize',40,'Linewidth',2)
plot([1 3 5],avAmp([16 21 20]),'.-','MarkerSize',40,'Linewidth',2)
plot(3,avAmp(18),'.-','MarkerSize',40,'Linewidth',2)
xlim([0.5 5.5]);
xticks([1:1:5]);
xticklabels({'12.5','25','50','75','87.5'})
legend('10','5','2','5mov','2mov','10mov')
xlabel('Duty Cycle ratio')
ylabel('amplitude')
saveas(gcf,'ampDC.jpg')

figure;hold on;
plot(static(1,:),'-^','MarkerSize',15,'LineWidth',3);plot(static(2,:),'-^','MarkerSize',15,'LineWidth',3);plot(static(3,:),'-^','MarkerSize',15,'LineWidth',3)
plot([2 3 4],moving([5 8 11]),'.--','MarkerSize',50,'Linewidth',3)
plot([1 3 5],moving([3 9 15]),'.:','MarkerSize',50,'Linewidth',1)
plot(3,moving(7),'.--','MarkerSize',50,'Linewidth',3)
xlim([0.5 5.5]);
xticks([1:1:5]);
ylim([0 3])
xticklabels({'12.5','25','50','75','87.5'})
legend('10','5','2','5mov','2mov','10mov','Location','BestOutside')
xlabel('Duty Cycle')
ylabel('ratings')
set(gca,'FontSize',20)
saveas(gcf,'ratingsDC.png')

%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Electrode based plots
%%%%%%%%%%%%%%%%%%

% plot mean amplitude for the 5 duty cycle at the 5 harmonics
% in 3 graphs separately for the 3 frequencies

pickElec = [23 15 28]; % 23 Oz, 9, B6=38, 16=best SNR
for ee=1:length(pickElec)
for cond=1:15 % exclude motion conditions
    for ff = 1 : 5
%         avData(cond).freq(avData(cond).i1f1*ff)
        meanAmp(cond,ff)=avData(cond).amp(pickElec(ee),(avData(cond).i1f1-1)*ff+1);
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
saveas(gcf,[num2str(pickElec(ee)) ' Amplitude.jpg'])

limCAx = [2.5 2.5 1 1];
figure;
for ff=1:4
subplot(2,2,ff)
imagesc(transpose([meanAmp(1:5,ff),meanAmp(6:10,ff),meanAmp(11:15,ff)]));
colorbar;caxis([0 limCAx(ff)]); 
xticks(1:5)
xticklabels({'12.5% DC','25% DC','50% DC','75% DC','87.5% DC'})
yticks(1:3)
yticklabels({'10.6Hz','5.3Hz','2.7Hz'})
title([num2str(pickElec(ee)) ' amplitude ' num2str(ff) 'f1'])
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [200 200 1000 800])
saveas(gcf,[num2str(pickElec(ee)) 'AmplitudeConfMat.jpg'])

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
saveas(gcf,[num2str(pickElec(ee)) ' TimeON.jpg'])


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
saveas(gcf,[num2str(pickElec(ee)) 'TimeOff.jpg'])

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
saveas(gcf,[num2str(pickElec(ee)) 'DutyCycle.jpg'])

% figure;line(8/85*1000 - (1/85*vect*1000),meanAmp(1:5,1),'Color','b','Marker','o')
% hold on; line(16/85*1000 - (1/85*2*vect*1000),meanAmp(6:10,1),'Color','r','Marker','o')
% hold on; line(32/85*1000 - (1/85*4*vect*1000),meanAmp(11:15,1),'Color','g','Marker','o')
% line(8/85*1000 - (1/85*vect*1000),meanAmp(1:5,2),'LineStyle','--','Color','b','Marker','*')
% hold on; line(16/85*1000 - (1/85*2*vect*1000),meanAmp(6:10,2),'LineStyle','--','Color','r','Marker','*')
% hold on; line(32/85*1000 - (1/85*4*vect*1000),meanAmp(11:15,2),'LineStyle','--','Color','g','Marker','*')
% xlabel('STA')

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% waveforms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pickElec = 23; % 23 Oz, 9, B6=38
% color = ['c','b','g','r','m'];
color = [[0 0 1];[0 0.8 1];[0 0.8 0];[0.8 0 0.8];[1 0 0]];

figure;hold on;
subplot(3,1,1);hold on;
for cond = 1:5
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'Color',color(cond,:),'LineWidth',3)
end
ylim([-2 2]);
legend('12','25','50','75','87','Location','BestOutside')
set(gca,'FontSize',20)
title('10Hz')
subplot(3,1,2);hold on;
for cond = 6:10
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'Color',color(cond-5,:),'LineWidth',3)
end
ylim([-4 5]);
legend('12','25','50','75','87','Location','BestOutside')
set(gca,'FontSize',20)
title('5Hz')
subplot(3,1,3);hold on;
for cond = 11:15
    plot(avData(cond).time,avData(cond).wave(pickElec,:),'Color',color(cond-10,:),'LineWidth',3)
end
ylim([-4 5]);
legend('12','25','50','75','87','Location','BestOutside')
set(gca,'FontSize',20)
title('2Hz')
saveas(gcf,'waveformOz.jpg')


figure;hold on;
subplot(3,1,1);hold on;
color = [0 0.8 0];
cc = [18];
for cond = 1:length(cc)
    plot(avData(cc(cond)).time,avData(cc(cond)).wave(pickElec,:),'Color',color(cond,:),'LineWidth',3)
end
ylim([-2 2]);
legend('50','Location','BestOutside')
set(gca,'FontSize',20)
title('5Hz')
color = [[0 0.8 1];[0 0.8 0];[0.8 0 0.8]];
cc = [17 22 19];
subplot(3,1,2);hold on;
for cond = 1:length(cc)
    plot(avData(cc(cond)).time,avData(cc(cond)).wave(pickElec,:),'Color',color(cond,:),'LineWidth',3)
end
ylim([-4 5]);
legend('25','50','75','Location','BestOutside')
set(gca,'FontSize',20)
title('2Hz')
color = [[0 0 1];[0 0.8 0];[1 0 0]];
cc = [16 21 20];
subplot(3,1,3);hold on;
for cond = 1:length(cc)
    plot(avData(cc(cond)).time,avData(cc(cond)).wave(pickElec,:),'Color',color(cond,:),'LineWidth',3)
end
ylim([-4 5]);
legend('12','50','87','Location','BestOutside')
title('1Hz')
set(gca,'FontSize',20)
saveas(gcf,'waveformMotion.jpg')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channel = 23;
figure;
for cond=1:15
    subplot(3,5,cond); hold on;
    bar(squeeze(avData(cond).freq((avData(cond).freq<15.5 & avData(1).freq>0))),squeeze(avData(cond).amp(channel,2:36)))
    if cond<6
        ylim([0 1.2]);
    else
        ylim([0 2.5]);
    end
end
set(gcf, 'Position', [0 500 1000 500])
saveas(gcf,'freqStaticOz.jpg')


position = [11 7 3 9 15 13 8];
figure;
for cond=16:22
    subplot(3,5,position(cond-15)); hold on;
    bar(squeeze(avData(cond).freq((avData(cond).freq<15.5 & avData(1).freq>0))),squeeze(avData(cond).amp(channel,2:36)))
    if cond==18
        ylim([0 1.2]);
    else
        ylim([0 2.5]);
    end
end
set(gcf, 'Position', [0 500 1000 500])
saveas(gcf,'freqMotionOz.jpg')
