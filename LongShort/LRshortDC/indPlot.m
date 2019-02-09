
% check ind sbj

dataDir = '/Users/marleneponcet/Documents/data/LRshortDC/V2/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

elec=23;
for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    
    figure; 
    subplot(2,2,1);hold on;
    plot(Axx(1).wave(elec,:))
    plot(Axx(2).wave(elec,:))
    plot(Axx(5).wave(elec,:))
    legend('AM','left','halfLeft')
    title('long range')
    subplot(2,2,2);hold on;
    bar(Axx(1).freq(2:30),Axx(1).amp(23,2:30))
    title('AM freq amplitude')
    
    subplot(2,2,3);hold on;
    plot(Axx(7).wave(elec,:))
    plot(Axx(8).wave(elec,:))
    plot(Axx(11).wave(elec,:))
    legend('AM','left','halfLeft')
    title('short range')
    subplot(2,2,4);hold on;
    bar(Axx(7).freq(2:30),Axx(7).amp(23,2:30))
    title('AM freq amplitude')    
    
    saveas(gcf,['sbj' num2str(ff)],'png')
end

for cond=1:size(sbjDiff,2)
    figure;hold on;
for ss=7%size(sbjDiff,1)
    plot(sbjDiff(ss,cond).data.wave(elec,:))
end
end


avPredictions = averageSbj(sbj([1:6 8:14],:));
elec=23; saveplot=0;
%%%%%% plot predictions
pickCond = [2:5 1; 7:10 6];
for ee=1:length(elec)
    for ss=1:2 % short and long range
        figure;hold on;
        for cc=1:length(pickCond)
            plot(avPredictions(pickCond(ss,cc)).time,avPredictions(pickCond(ss,cc)).wave(elec(ee),:),'LineWidth',2);
        end
        ylim([-4 6])
        plot([0 400],[0 0],'--k')
        xlabel('time (ms)')
        ylabel('amplitude')
        title(['electrode' num2str(elec(ee))])
        legend('linear','spatial','temp','spatial+temp','AM')
        if saveplot
            saveas(gcf,[labelIn{dd} 'predict' titolo{ss} num2str(elec(ee)) '.pdf'],'pdf')
            saveas(gcf,[labelIn{dd} 'predict' titolo{ss} num2str(elec(ee)) '.fig'],'fig')
        end
    end
end

%%% differences
avDiff = averageSbj(sbjDiff([1:6 8:14],:));
for ee=1:length(elec)
    for ss=1:2 % short and long range
        figure;hold on;
        for cc=1:4
            plot(avDiff(cc+4*(ss-1)).time,avDiff(cc+4*(ss-1)).wave(elec(ee),:),'LineWidth',2);
        end
        ylim([-3 6])
        plot([0 400],[0 0],'--k')
        xlabel('time (ms)')
        ylabel('amplitude difference')
        title(['Diff for electrode' num2str(elec(ee))])
        legend('linear','spatial','temp','spatial+temp')
        if saveplot
            saveas(gcf,[labelIn{dd} 'predict' titolo{ss} num2str(elec(ee)) '.pdf'],'pdf')
            saveas(gcf,[labelIn{dd} 'predict' titolo{ss} num2str(elec(ee)) '.fig'],'fig')
        end
    end
end
 
    figure; hold on;
for cc=6:10    
plot(avPredictions(cc).time,avPredictions(cc).wave(elec(ee),:),'LineWidth',2);
end
legend('AM','linear','spatial','temp','spatial+temp')