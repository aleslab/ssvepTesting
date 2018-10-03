close all;
clearvars;
% ATTENTION!!!!!!!!!!!
% pour exp 1, le cyle est inverse par rapport a exp 2 (rien puis single
% flash ou rien puis simult flashes). Change time axis

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% addpath C:\Users\Marlene\Documents\git\fieldtrip-aleslab-fork
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\svndlCopy
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\biosemiUpdated
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFuntions
% ft_defaults
% dataIn = {'C:\Users\Marlene\Documents\git\dataLR\LRlongDC\', 'C:\Users\Marlene\Documents\git\dataLR\LRshortDC\V2\'};


dataIn = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};
labelIn = {'75DC','25DC'};
colDiff = {'r','g','b','c'};
elec = [23]; % 23=Oz, 19=Pz, 3=CPz, A7/B4-36=P3/P4, FCz=87, C3=115, C4=54
saveplot = 1;

for dd=1:length(dataIn) % 2 experiments
    clear cfg avPredictions sbj diff avInteractions interaction;
    %%% load the sbj predictions
    load([dataIn{dd} 'sbjprediction.mat'])
    load([dataIn{dd} 'NLinteraction.mat'])
    
    % do the average
    if dd==1
    avPredictions = averageSbj(sbj);
    avInteractions = averageSbj(interaction);
    else
    avPredictions = averageSbj(sbj([1:6 8:end],:));
    avInteractions = averageSbj(interaction([1:6 8:end],:));
    end
    
    %%%%
    avPredictions(1).condLabel = 'originalmotion';
    avPredictions(2).condLabel = 'linear';
    avPredictions(3).condLabel = 'spatial';
    avPredictions(4).condLabel = 'temp';
    avPredictions(5).condLabel = 'spatiotemp';
    % avPredictions(6).condLabel = 'originalmotion';
    % avPredictions(7).condLabel = 'SR linearPred';
    % avPredictions(8).condLabel = 'SR spatialPred';
    % avPredictions(9).condLabel = 'SR tempPred';
    % avPredictions(10).condLabel = 'SR spatiotempPred';
    % avPredictions(11).condLabel = 'LR STnl';
    % avPredictions(12).condLabel = 'SR STnl';

    
    % create a new field "filteredWave" after low-pass filter
    % warning about filter is NORMAL
    for condIdx=1:length(avPredictions)
        filtIdx = determineFilterIndices( 'low49', avPredictions(condIdx).freq, avPredictions(condIdx).i1f1 );
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(avPredictions(condIdx).amp));
        filtMat(:,filtIdx) = true;
        
        %Combine the filter and sig vaules with a logical AND.
        % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        avPredictions(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
        cfg.activeFreq =  avPredictions(condIdx).activeFreq;
        
        [ avPredictions(condIdx).filteredWave ] = filterSteadyState( cfg, avPredictions(condIdx) );
    end
    for condIdx=1:length(avInteractions)
        filtIdx = determineFilterIndices( 'low49', avInteractions(condIdx).freq, avInteractions(condIdx).i1f1 );
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(avInteractions(condIdx).amp));
        filtMat(:,filtIdx) = true;
        
        %Combine the filter and sig vaules with a logical AND.
        % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        avInteractions(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
        cfg.activeFreq =  avInteractions(condIdx).activeFreq;
        
        [ avInteractions(condIdx).filteredWave ] = filterSteadyState( cfg, avInteractions(condIdx) );
    end
    
    
    %%% Correction shift cycle only for E1
    if dd==1
    for cond = 1:length(avPredictions)
        avPredictions(cond).filteredWave = circshift(avPredictions(cond).filteredWave,[0 length(avPredictions(1).filteredWave)/2]);
    end    
    for cond =1:length(avInteractions)
       avInteractions(cond).filteredWave = circshift(avInteractions(cond).filteredWave,[0 length(avInteractions(1).filteredWave)/2]);
    end
    end   
    
    
    %%%%% plot long and short range
    for ee=1:length(elec)
        figure;hold on;
        plot(avPredictions(1).time,avPredictions(1).filteredWave(elec(ee),:),'LineWidth',2);
        plot(avPredictions(6).time,avPredictions(6).filteredWave(elec(ee),:),'LineWidth',2);
        ylim([-3 3])
        plot([0 400],[0 0],'--k')
        xlabel('time (ms)')
        ylabel('amplitude')
        title(['electrode' num2str(elec(ee))])
        legend('LR','SR')
        if saveplot
            saveas(gcf,[labelIn{dd} '_LRSR_' num2str(elec(ee)) '.pdf'],'pdf')
            saveas(gcf,[labelIn{dd} '_LRSR_' num2str(elec(ee)) '.fig'],'fig')
        end
    end
    
    
    %%%%% plot interactions
    titolo = {'LR','SR'};
    for ee=1:length(elec)
        for ss=1:2 % short and long range
            figure;hold on;
            for cc=1:3
                plot(avInteractions(cc+3*(ss-1)).time,avInteractions(cc+3*(ss-1)).filteredWave(elec(ee),:),colDiff{cc},'LineWidth',2);
            end
            plot([0 400],[0 0],'--k')
            ylim([-3 3])
            xlabel('time (ms)')
            ylabel('amplitude')
            title(['electrode' num2str(elec(ee)) titolo{ss}])
            legend('spatial','templeft','tempright')
            if saveplot
                saveas(gcf,[labelIn{dd} 'interaction' titolo{ss} num2str(elec(ee)) '.pdf'],'pdf')
                saveas(gcf,[labelIn{dd} 'interaction' titolo{ss} num2str(elec(ee)) '.fig'],'fig')
            end
        end
        
    end
    
    
    
    %%%%%% plot predictions
    pickCond = [2:5 1 ; 7:10 6];
    for ee=1:length(elec)
        for ss=1:2 % short and long range
        figure;hold on;
        for cc=1:length(pickCond)
            plot(avPredictions(pickCond(ss,cc)).time,avPredictions(pickCond(ss,cc)).filteredWave(elec(ee),:),'LineWidth',2);
        end
        ylim([-3 6])
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
    
    
    
    %%%%%% compute differences
    pickCond = [2:5 ; 7:10];
    for cc=1:size(pickCond,1)
        condCompare = pickCond(cc,:);
        fixCond = condCompare(1) - 1;
        for cond=1:length(condCompare)
            diff(cc,cond) = computeDiff(avPredictions(fixCond),avPredictions(condCompare(cond)));
        end
    end
    % create a new field "filteredWave" after low-pass filter
    % warning about filter is NORMAL
    for cc=1:2 % LR - SR
        for condIdx=1:length(diff)
            filtIdx = determineFilterIndices( 'low49', diff(cc,condIdx).freq, diff(cc,condIdx).i1f1 );
            
            %Create a logical matrix selecting frequency components.
            filtMat = false(size(diff(cc,condIdx).amp));
            filtMat(:,filtIdx) = true;
            
            %Combine the filter and sig vaules with a logical AND.
            % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
            diff(cc,condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
            cfg.activeFreq =  diff(cc,condIdx).activeFreq;
            
            [ diff(cc,condIdx).filteredWave ] = filterSteadyState( cfg, diff(cc,condIdx) );
        end
    end
    
    % plot the wave differences
    for ee=1:length(elec)
        for cc=1:2 % LR - SR
            figure; hold on;
            for cond=1:4
                plot(diff(cc,cond).time,diff(cc,cond).filteredWave(elec(ee),:),'LineWidth',2);
            end
            plot([0 400],[0 0],'--k')
            ylim([-3 3])
            title([labelIn{dd} 'diff' titolo{cc}])
            xlabel('time (ms)')
            ylabel('difference')
            legend(avPredictions(2:5).condLabel)
            if saveplot
                saveas(gcf,[labelIn{dd} 'diff' titolo{cc} '-e' num2str(elec(ee)) '.fig'],'fig')
                saveas(gcf,[labelIn{dd} 'diff' titolo{cc} '-e' num2str(elec(ee)) '.pdf'],'pdf')
            end
        end
        
    end
    
    
    % plot the freq (check if it is all noise..)
    maxFreq = 15;
    for ee=1:length(elec)
        for cc=1:2 % LR - SR
            figure; hold on;
            indexF = diff(cc,1).freq(diff(cc,1).freq<maxFreq & diff(cc,1).freq>0);
            bar(indexF,squeeze(diff(cc,1).amp(elec(ee),2:length(indexF)+1)),0.4,'FaceColor','r','LineStyle','none')
            bar(indexF,squeeze(diff(cc,4).amp(elec(ee),2:length(indexF)+1)),0.6,'FaceColor','none','LineWidth',1.5,'EdgeColor','c')
            ylim([0 1.1])
            title([labelIn{dd} 'freq' titolo{cc} '-e' num2str(elec(ee))])
            if saveplot
                saveas(gcf,[labelIn{dd} 'freq' titolo{cc} '-e' num2str(elec(ee)) '.fig'],'fig')
                saveas(gcf,[labelIn{dd} 'freq' titolo{cc} '-e' num2str(elec(ee)) '.pdf'],'pdf')
            end
        end
        
    end
    
    % plot topo freq diff
    %     harm = avPredictions(1).freq(avPredictions(1).i1f1).*[1:6];
    harm = avPredictions(1).freq(avPredictions(1).i1f1).*[2 4];
    for hh=1:length(harm)
        harmIndex(hh) = find(avPredictions(1).freq==harm(hh));
    end
    for lr=1:2
    for cond=[1 4]
    figure; hold on;
    for hh=1:length(harmIndex)
        subplot(1,2,hh)
        plotTopo(diff(lr,cond).amp(:,harmIndex(hh)),cfg.layout);
        caxis([0 1])
        colorbar
        title([titolo{lr} num2str(harm(hh)) 'Hz'])
    end
    if saveplot
        saveas(gcf,[labelIn{dd} 'topofreq Cond' num2str(cond) titolo{lr} num2str(harm(hh)) 'Hz' '.fig'],'fig')
        saveas(gcf,[labelIn{dd} 'topofreq Cond' num2str(cond) titolo{lr} num2str(harm(hh)) 'Hz' '.pdf'],'pdf')
    end
    end
    end
end


% figure;hold on;
% for ss=1:7
%     test(ss,:)=(sbj(ss,6).data.wave(23,:)- sbj(ss,5).data.wave(23,:));
%     plot((sbj(ss,1).data.wave(23,:)- sbj(ss,5).data.wave(23,:)))
% end
% plot(mean(test(:,:),1),'r','LineWidth',3)
% plot(diff(1,2).wave(23,:),'g','LineWidth',3)


