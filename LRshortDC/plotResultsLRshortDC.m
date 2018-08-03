
clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/data/LRshortDC/V2/Axx/';
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
subplot(2,2,ss); hold on;
plotTopo(snr(:,ss),cfg.layout);
colorbar; 
% caxis([1 4.5]);
end
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 0 2000 2500])
saveas(gcf,'snrSbj.jpg')


%% average
keepSbj = 1:4;

avData = averageAxx(dataSbj(:,keepSbj));

    aa = {'Long','Short'}; bb=[0 6];
    for tt=1:2
    avData(1+bb(tt)).condLabel = [aa{tt} 'Motion'];
    avData(2+bb(tt)).condLabel = [aa{tt} 'Left'];
    avData(3+bb(tt)).condLabel = [aa{tt} 'Right'];
    avData(4+bb(tt)).condLabel = [aa{tt} 'Simult'];
    avData(5+bb(tt)).condLabel = [aa{tt} 'HalfLeft'];
    avData(6+bb(tt)).condLabel = [aa{tt} 'HalfRight'];
    end
    
% compute some stats
% avData = gpStatT2_DC(dataSbj,avData);

interactiveSteadyStatePlot2(cfg,avData) 



%% 
pickElec = [23 8 37];
color = [[0 0 1];[0 0.8 1];[0 0.8 0];[0.8 0 0.8];[1 0 0];[0 0 0]];
for pp=1:length(pickElec)
figure;hold on;
for ss=1:size(dataSbj,2)
subplot(2,4,1+2*(ss-1)); hold on;
for cond=1:6
    plot(dataSbj(cond,ss).wave(pickElec(pp),:),'Color',color(cond,:),'LineWidth',2)
end
legend('LR','L','R','simult','HL','HR','Location','Best')
subplot(2,4,2+2*(ss-1)); hold on;
for cond=7:12
    plot(dataSbj(cond,ss).wave(pickElec(pp),:),'Color',color(cond-6,:),'LineWidth',2)
end
legend('SR','L','R','simult','HL','HR','Location','Best')
end
end

for pp=1:length(pickElec)
figure;hold on;
subplot(2,1,1); hold on;
for cond=1:6
    plot(avData(cond).time,avData(cond).wave(pickElec(pp),:),'Color',color(cond,:),'LineWidth',2)
end
legend('LR','L','R','simult','HL','HR','Location','BestOutside')
subplot(2,1,2); hold on;
for cond=7:12
    plot(avData(cond).time,avData(cond).wave(pickElec(pp),:),'Color',color(cond-6,:),'LineWidth',2)
end
legend('SR','L','R','simult','HL','HR','Location','BestOutside')
saveas(gcf,['Figure' filesep num2str(pickElec(pp)) 'WavePerCond.png'],'png')
end



%%%
load('sbjprediction.mat')
elec = [23 8 37];
for pick=1:length(elec)
figure;
for ss=1:4
    subplot(2,2,ss);hold on;
    for cond=1:5
        plot(sbj(ss,cond).data.wave(elec(pick),:),'LineWidth',3)
    end
    legend('LR','linear','lin+spatial','lin+temp','lin+spat+temp','Location','Best')
end
end
for pick=1:length(elec)
figure;
for ss=1:4
    subplot(2,2,ss);hold on;
    for cond=6:10
        plot(sbj(ss,cond).data.wave(elec(pick),:),'LineWidth',3)
    end
    legend('SR','linear','lin+spatial','lin+temp','lin+spat+temp','Location','Best')
end
end

gpAve = averageSbj(sbj);
for pick=1:length(elec)
    figure;
    subplot(2,1,1);hold on;
    for cond=1:5
        plot(gpAve(cond).time,gpAve(1,cond).wave(elec(pick),:),'LineWidth',3)
    end
    legend('LR','linear','lin+spatial','lin+temp','lin+spat+temp','Location','BestOutside')
    subplot(2,1,2);hold on;
    for cond=6:10
        plot(gpAve(cond).time,gpAve(1,cond).wave(elec(pick),:),'LineWidth',3)
    end
    legend('SR','linear','lin+spatial','lin+temp','lin+spat+temp','Location','BestOutside')
    saveas(gcf,['Figure' filesep num2str(elec(pick)) 'WavePred.png'],'png')
end


load('NLinteraction.mat')
for pick=1:length(elec)
figure;
for ss=1:4
    subplot(2,2,ss);hold on;
    for cond=1:8
        plot(interaction(ss,cond).wave(elec(pick),:),'LineWidth',3)
    end
    legend('LRspatial','LRl','LRr','SRspatial','SRl','SRr','LRst','SRst','Location','Best')
end
end
gpInteraction = averageAxx(permute(interaction,[2 1]));
for pick=1:length(elec)
figure;
    subplot(2,1,1);hold on;
    for cond=[1 2 3 7]
        plot(gpInteraction(cond).time,gpInteraction(cond).wave(elec(pick),:),'LineWidth',3)
    end
    legend('LRspatial','LRl','LRr','LRst','Location','BestOutside')
    subplot(2,1,2);hold on;
    for cond=[4 5 6 8]
        plot(gpInteraction(cond).time,gpInteraction(cond).wave(elec(pick),:),'LineWidth',3)
    end
    legend('SRspatial','SRl','SRr','SRst','Location','BestOutside')
    saveas(gcf,['Figure' filesep num2str(elec(pick)) 'Interaction.png'],'png')
end






%% freq plots
condName={'LR','L','R','simult','HL','HR'};
for pp=1:length(pickElec)
figure;hold on;
for cond=1:6
    subplot(3,2,cond); hold on;
    bar(avData(cond).freq(1:37),avData(cond).amp(pickElec(pp),1:37))
    title(condName(cond));
end
figure;hold on;
for cond=7:12
    subplot(3,2,cond-6); hold on;
    bar(avData(cond).freq(1:37),avData(cond).amp(pickElec(pp),1:37))
end
end
