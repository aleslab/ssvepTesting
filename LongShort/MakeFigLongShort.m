% Make figures
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% see pred.labels for details
% first 12 are 1st exp, from 12 is the 2nd exp
load('results.mat')
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

% find harmonics for SSVEPs
harmIndex = [];
for pp = [1 11] % go through the 2 exp
    harm = pred(pp).freq(pred(pp).i1f1).*[1 2 3 4];
    for hh=1:length(harm)
        harmI(hh) = find(pred(pp).freq==harm(hh));
    end
    harmIndex = [harmIndex; harmI];
end
% index 11 & 21 (E1) 13 & 25 (E2)

elec = [23 9 38];

%%%%% plot long and short range
for ee=1:2
    figure;hold on;
    ff=0;
    for pp = 1:5:10
        for hh = 1:length(harmIndex) % look at the 2 mains harmonics
            ff=ff+1;
            subplot(2,4,ff)
            plotTopo(pred(pp+(ee-1)*10).amp(:,harmIndex(ee,hh)),cfg.layout);
            %         caxis([0 2.5]);
            colorbar;
            title(['E' num2str(ee) '-1F' num2str(hh)])
        end
    end
end


% figure;hold on;
% ff=0;
% for pp = 1:5:20 % go through the 2 exp
%     if pp< 11
%         ee=1;
%     else
%         ee=2;
%     end
%     ff=ff+1;
%     subplot(2,2,ff)
%     plotTopo(pred(pp).amp(:,harmIndex(ee,2)),cfg.layout);
%     caxis([0 2.5]);
% %     colorbar;
%     title(['E' num2str(ee) '-1F2'])
% end
% saveas(gcf,'topoAM.pdf','pdf')


% for ee=1:2
% figure;hold on;
% for chan=1:length(elec)
%     subplot(1,3,chan); hold on;
%     for pp=1:5:10
%         plot(pred(pp+10*(ee-1)).time,pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
%     end
%     title(elec(chan))
%     legend('LR','SR')
% end
% end

% %%%%%%% plot topo averaged over 400 ms (1 cycle)
% ppp=0;
% figure;hold on;
% for ee=1:2 % experiment 1 and 2
% for pp=1:5:10 % AM long and short
%     ppp = ppp+1;
%     subplot(2,2,ppp); hold on;
%     plotTopo(mean(abs(pred(pp+(ee-1)*10).filteredWave),2),cfg.layout);
%     colorbar;
% %     caxis([0 1.5])
% end
% end


for chan=1:length(elec)
    subplot(1,3,chan); hold on;
    for pp=1:5:10
        plot(pred(pp+10*(ee-1)).time,pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
    end
    title(elec(chan))
    legend('LR','SR')
end


%%%%%%% plot predictions
for ee=1:2
figure('Renderer', 'painters', 'Position', [10 10 1200 500])
for chan=1:length(elec)
    subplot(2,3,chan); hold on;
    for pp=1:5
        plot(pred(pp+10*(ee-1)).time,pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
    end
    title(['LR' num2str(elec(chan))])
    legend('AM','linear','spatial','temp','ST','Location','Best')
    subplot(2,3,chan+3); hold on;
    for pp=6:10
        plot(pred(pp+10*(ee-1)).time,pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
    end
    title(['SR' num2str(elec(chan))])
    legend('AM','linear','spatial','temp','ST','Location','Best')
end
saveas(gcf,['E' num2str(ee) 'predictions'],'png')
end

%%%%%%% plot differences
for ee=1:2
figure('Renderer', 'painters', 'Position', [10 10 1200 500])
for chan=1:length(elec)
    subplot(2,3,chan); hold on;
    for pp=1:5
        plot(pred(pp+10*(ee-1)).time,pred(1+10*(ee-1)).filteredWave(elec(chan),:)-pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
    end
    title(['LR' num2str(elec(chan))])
    legend('AM','linear','spatial','temp','ST','Location','Best')
    subplot(2,3,chan+3); hold on;
    for pp=6:10
        plot(pred(pp+10*(ee-1)).time,pred(6+10*(ee-1)).filteredWave(elec(chan),:)-pred(pp+10*(ee-1)).filteredWave(elec(chan),:),'LineWidth',2);
    end
    title(['SR' num2str(elec(chan))])
    legend('AM','linear','spatial','temp','ST','Location','Best')
end
saveas(gcf,['E' num2str(ee) 'differences'],'png')
end

    
  
%%%%%%%%% DIFFERENCES    
load('diff.mat')
elec = [23 9 38];


% RMS root mean square
% diff(exp condition, factor tested)
% exp condition: LRE1 SRE1 LRE2 SRE2
% factor tested: linear spatial temporal ST
for ec=1:4
    for fact=1:4
%          rmsElec = rms(diff(ec, fact).wave,2); % rms per electrode
%          rootMeanS(ec,fact,:) = rmsElec
         rootMeanS(ec,fact) = rms(diff(ec, fact).wave(:));
    end
end
figure; bar(rootMeanS)
legend({'linear' 'spatial' 'temporal' 'ST'},'location','Best')
xticklabels({'LR1' 'SR1' 'LR2' 'SR2'})
ylabel('mean RMS')
saveas(gcf,'meanRMS','png')


%%%%%%% RMS TOPO
for ec=1:4
    for fact=1:4
        for trode = 1:size(diff(ec, fact).wave(:,:),1)
            rootMeanElec(trode,ec,fact) = rms(diff(ec, fact).wave(trode,:));
        end
    end
end

for ec=1:4
    figure;hold on;
    for fact=1:4
        subplot(2,2,fact); hold on;
        plotTopo(rootMeanElec(:,ec,fact),cfg.layout);
        colorbar;
    end
end

%%%%%%% RMS per electrode
for trode = 1:length(elec)
    figure; bar(squeeze(rootMeanElec(elec(trode),:,:)))
end





% test = [2 2; 1 1; 3 3];
% rms(test,2)
% sqrt(mean(test .* test,2))
% rms(rms(test,2))
% rms(test(:))
% 
% rms(diff(ec, fact).wave(2:4,2:3),2)
% sqrt(mean(diff(ec, fact).wave(2:4,2:3) .* diff(ec, fact).wave(2:4,2:3),2))



% for ee=1:2
%     figure('Renderer', 'painters', 'Position', [10 10 1200 500])
%     for chan=1:length(elec)
%         subplot(2,3,chan); hold on;
%         for cc=1:4
%             plot(diff(1+2*(ee-1),cc).time,diff(1+2*(ee-1),cc).filteredWave(elec(chan),:),'LineWidth',2);
%         end
%         title([num2str(elec(chan))])
%         legend('linear','spatial','temp','ST','Location','Best')
%         subplot(2,3,chan+3); hold on;
%         for cc=1:4
%             plot(diff(2+2*(ee-1),cc).time,diff(2+2*(ee-1),cc).filteredWave(elec(chan),:),'LineWidth',2);
%         end
%         title([num2str(elec(chan))])
%         legend('linear','spatial','temp','ST','Location','Best')
%     end
% saveas(gcf,['E' num2str(ee) 'differences2'],'png')
% end
% 
% 
% 
% %%%% plot SSVEP topo
% for ee=1:2
%     figure;hold on;
%     ff=0;
%     for pp = 1:5:10
%         for hh = 1:length(harmIndex) % look at the 2 mains harmonics
%             ff=ff+1;
%             subplot(2,4,ff)
%             plotTopo(diff(pp+(ee-1)*10).amp(:,harmIndex(ee,hh)),cfg.layout);
%             %         caxis([0 2.5]);
%             colorbar;
%             title(['E' num2str(ee) '-1F' num2str(hh)])
%         end
%     end
% end
% 
% % plot the freq (check if it is all noise..)
% maxFreq = 15;
% titolo = {'linear','spatial','temp','ST'};
% type = {'LR','SR','LR','SR'};
% exper = [1 1 2 2];
% for chan=1:length(elec)
%     figure('Renderer', 'painters', 'Position', [10 10 1200 500])
%     for cc=1:4
%         indexF = diff(cc,1).freq(diff(cc,1).freq<maxFreq & diff(cc,1).freq>0);
%         for inter = 1:4
%             subplot(4,4,inter+(cc-1)*4)
%             bar(indexF,squeeze(diff(cc,inter).amp(elec(chan),2:length(indexF)+1)),0.4,'FaceColor','r','LineStyle','none')
%             ylim([0 1.4])
%             title(['E' num2str(exper(cc)) type{cc} titolo{inter} 'Chan' num2str(elec(chan))])
%         end
%     end
%     saveas(gcf,['leftOverElec' num2str(elec(chan))],'png')
% end



