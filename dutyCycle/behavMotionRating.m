
clearvars
close all

dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/ratings/';

behavFiles = dir([dataDir 'DC*.mat']);

for ff=1:length(behavFiles)
    clear strMot stdTab avTab allR;
    load([dataDir behavFiles(ff).name])

    %remove any invalid trial from the data
    if length(experimentData) ~= 198 || ff==1 % this sbj has different nb of trials
        toRemove = [];
        for ll=1:length(experimentData)
            if experimentData(ll).validTrial == 0 || isempty(experimentData(ll).trialData)
                toRemove = [toRemove; ll];
            end
        end
    else
        toRemove = [];
    end
    keepIndexes = setdiff(1:length(experimentData), toRemove');
    tableData = struct2table(experimentData(keepIndexes));
    
    
    allCond = unique(tableData.condNumber);
    
    for cc=1:length(allCond)
        clear arrResp;
        indexCond = find(tableData.validTrial & tableData.condNumber == allCond(cc));
        for ic = 1: length(indexCond)
            arrResp(ic)  = tableData.trialData(indexCond(ic)).response;
        end
        allR(cc,:) = arrResp;
        stdTab(cc) = std(arrResp);
        avTab(cc) = mean(arrResp);
    end

    strMot = avTab;
    tabFlash = [strMot(1:5);strMot(6:10);strMot(11:15)];
    tabPmot = zeros(3,5);
    tabPmot(3,1)= strMot(16);
    tabPmot(2,2)= strMot(17);
    tabPmot(1,3)= strMot(18);
    tabPmot(2,4)= strMot(19);
    tabPmot(3,5)= strMot(20);
    tabPmot(3,3)= strMot(21);
    tabPmot(2,3)= strMot(22);
    tabStatic(:,:,ff) = tabFlash;
    tabMot(:,:,ff) = tabPmot;

end

% each ind
for ss=1:length(keepSBJ)
    figure;imagesc(tabStatic(:,:,ss))
    figure;imagesc(tabMot(:,:,ss))
end

figure;
subplot(1,2,1)
imagesc(mean(tabStatic,3),[0 3]);
subplot(1,2,2)
imagesc(mean(tabMot,3),[0 3]);colorbar;
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Position', [0 1000 1200 500])
saveas(gcf,'figures/ratingsDC.jpg')

static = mean(tabStatic,3); moving = mean(tabMot,3);
save('ratings10.mat','static','moving','tabStatic','tabMot')


%%%%%%%%%%% PLOT
dcVal = [12.5 25 50 75 87.5];
for fq=1:3
    for dc=1:5
        statSEM(fq,dc) = std(tabStatic(fq,dc,:))/sqrt(length(tabStatic));
        movSEM(fq,dc) = std(tabMot(fq,dc,:))/sqrt(length(tabMot));
    end
end
figure;hold on;
for fq=1:3
    errorbar(dcVal,static(fq,:),statSEM(1,:),['.-'  col{fq}],'MarkerSize',40,'LineWidth',2)
end
errorbar(dcVal(3),moving(7),movSEM(7),['^'  col{1}],'MarkerSize',15,'Linewidth',2)
errorbar([25 50 75],moving([5 8 11]),movSEM([5 8 11]),['^:'  col{2}],'MarkerSize',15,'Linewidth',2)
errorbar([12.5 50 87.5],moving([3 9 15]),movSEM([3 9 15]),['^:'  col{3}],'MarkerSize',15,'Linewidth',2)
xlim([0 100]);
xticks([0:12.5:100]);
ylim([0 3])
legend('10','5','2.5','Location','Best')
xlabel('Duty Cycle')
ylabel('motion ratings')
set(gca,'FontSize',15)
saveas(gcf,['figures' filesep 'ratingsDC.png'])

% avMot = mean(tabMot,3);
% valInd = find(avMot>0);
% listX = [12.5 25 50 50 50 75 87.5];
% listCol = [3 2 1 2 3 2 3];
% col = {'b','r','g'};
% figure; hold on;
% for freq=1:3
%     plot([12.5 25 50 75 87.5], mean(tabStatic(freq,:,:),3),['--o' col{freq}],'LineWidth',2)
% end
% for vv=1:length(valInd)
%     scatter(listX(vv), avMot(valInd(vv)),80, 'filled',col{listCol(vv)})
% end
% xlim([0 100])
% xticks([12.5 25 50 75 87.5])
% ylim([0 3])
% xlabel('duty cycle')
% legend('10','5','2.5','location','best')
% saveas(gcf,'figures/ratingPlot.png')
