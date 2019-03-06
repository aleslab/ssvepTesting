% compute RMS

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

% see pred.labels for details
% first 12 are 1st exp, from 12 is the 2nd exp
clearvars;
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};

for ee=1:2 % which experiment
    clear diffFilt rmsFilt rmseNoCoef rmsNoise rmsSignal rmseCoef rmseNoCoef testrmsNoise rmsTopo rmsTopoNoise
    load([dataPath{ee} 'sbjprediction.mat'])
    load(['rmseCoefE' num2str(ee) '.mat'])
    
    %%%%%%% RMS root mean square
    % calculate per participant using filteredwave
    for subj=1:length(sbj)
        for chan=1:sbj(1,1).data.nchan
            for ss=1:2 % long/short
                for fact =1:4
                    diffFilt(fact+4*(ss-1),subj,chan,:) = (sbj(subj,fact+1+5*(ss-1)).data.filteredWave(chan,:) - sbj(subj,1+5*(ss-1)).data.filteredWave(chan,:));
                    rmsFilt(fact+4*(ss-1),subj,chan) = rms(diffFilt(fact+4*(ss-1),subj,chan,:));
                end
            end
        end
    end

    
    %%%%%%%%%%%%%%%%%
    % look at rms per condition
    for ss=1:length(sbj)
        for numCond=1:10
            testrmsNoise(ss,numCond,:) = rms(sbj(ss,numCond).data.noiseWave(:)  );
            rmsSignal(ss,numCond,:) = rms(sbj(ss,numCond).data.filteredWave(:)  );
        end
    end
    figure;hold on;
    bar(mean(testrmsNoise))
    errorbar(mean(testrmsNoise), std(testrmsNoise),'k','Linestyle','none')
    title(['E' num2str(ee)])
    ylabel('rms noise per cond')
    saveas(gcf,['figures' filesep 'RMSnoiseAllCond' num2str(ee)],'png')
    
    figure;hold on;
    bar(mean(rmsSignal))
    errorbar(mean(rmsSignal), std(rmsSignal),'k','Linestyle','none')
    title(['E' num2str(ee)])
    ylabel('rms Signal per cond')
    
    %%%%%%%
    % get all the rms we need across electrodes
    rmseNoCoef(:,:) = rms(rmsFilt(:,:,:),3)';
    rmsNoise = rms(testrmsNoise(:,:,:),3);
    save(['allRMSe' num2str(ee) '.mat'],'rmsSignal','rmseNoCoef','rmsNoise','rmseCoef')
    
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Topography rmsTopo AM, rms Lin, rmse Lin, rmse Temp
    % for this I need rms per electrode
    
    for ss=1:length(sbj)
        for chan = 1 : size(sbj(ss,numCond).data.filteredWave,1)
            rmsTopo(ss,1,chan) = rms(sbj(ss,1).data.filteredWave(chan,:) );
            rmsTopo(ss,2,chan) = rms(sbj(ss,2).data.filteredWave(chan,:) );
            rmsTopo(ss,3,chan) = rms(sbj(ss,2).data.filteredWave(chan,:) - sbj(ss,1).data.filteredWave(chan,:) );
            rmsTopo(ss,4,chan) = rms(sbj(ss,4).data.filteredWave(chan,:) - sbj(ss,1).data.filteredWave(chan,:) );
            
            rmsTopo(ss,5,chan) = rms(sbj(ss,6).data.filteredWave(chan,:) );
            rmsTopo(ss,6,chan) = rms(sbj(ss,7).data.filteredWave(chan,:) );
            rmsTopo(ss,7,chan) = rms(sbj(ss,7).data.filteredWave(chan,:) - sbj(ss,6).data.filteredWave(chan,:) );
            rmsTopo(ss,8,chan) = rms(sbj(ss,9).data.filteredWave(chan,:) - sbj(ss,6).data.filteredWave(chan,:) );
        end
    end
    
    % and I will need to divide by the noise!!! 
    for ss=1:length(sbj)
        for chan = 1 : size(sbj(ss,numCond).data.filteredWave,1)
            cc=0;
            for numCond=[1:4 6:9] % do not include pred s+t
                cc=cc+1;
                rmsTopoNoise(ss,cc,chan) = rms(sbj(ss,numCond).data.noiseWave(chan,:)  );
            end
        end
    end
    save(['chanRMS' num2str(ee) '.mat'],'rmsTopo','rmsTopoNoise')    
    
end





% figure;
% for ee=1:2
%     load(['allRMSe' num2str(ee)]);
%     subplot(2,2,1+2*(ee-1));hold on;
%     bar([mean(rmsSignal(:,1:2)) mean(rmsNoise(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))])
%     errorbar([mean(rmsSignal(:,1:2)) mean(rmsNoise(:,1:2)) mean(rmseNoCoef(:,1:4)) mean(rmseCoef(:,1:3))], ...
%         [std(rmsSignal(:,1:2)) std(rmsNoise(:,1:2)) std(rmseNoCoef(:,1:4)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
%     title(['E' num2str(ee) 'LR'])
%     subplot(2,2,2+2*(ee-1));hold on;
%     bar([mean(rmsSignal(:,6:7)) mean(rmsNoise(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))])
%     errorbar([mean(rmsSignal(:,6:7)) mean(rmsNoise(:,6:7)) mean(rmseNoCoef(:,5:8)) mean(rmseCoef(:,4:6))], ...
%         [std(rmsSignal(:,6:7)) std(rmsNoise(:,6:7)) std(rmseNoCoef(:,5:8)) std(rmseCoef(:,1:3))],'k','Linestyle','none')
%     title(['E' num2str(ee) 'SR'])
% end
% saveas(gcf,['figures' filesep 'RMSallCond' num2str(ee)],'png')


