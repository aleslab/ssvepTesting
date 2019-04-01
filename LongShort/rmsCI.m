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
    
    load([dataPath{ee} 'sbjprediction.mat'])
    
    % conditions in sbj: AM, linear, spatial, temp, S+T, AM short, linear,
    % spatial, temp, S+T, non-linear ST long, non-linear ST short
    
    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% over all time and electrodes
    clear rmsNoise rmsSignal rmse nrms nrmse nrmseCoef
    % RMS over all time and electrodes
    for ss=1:length(sbj)
        for numCond=1:10 
            rmsNoise(ss,numCond) = rms(sbj(ss,numCond).data.noiseWave(:)  );
            rmsSignal(ss,numCond) = rms(sbj(ss,numCond).data.filteredWave(:)  );
        end
    end
    % RMSE over all time and electrodes
    % lin, spa, temp, s+t
    for ss=1:length(sbj)
        for numCond=2:5 
            rmse(ss,numCond-1) = rms(sbj(ss,1).data.filteredWave(:) -  sbj(ss,numCond).data.filteredWave(:)); % long-range
            rmse(ss,numCond+3) = rms(sbj(ss,6).data.filteredWave(:) -  sbj(ss,numCond+5).data.filteredWave(:)); % short-range
        end
    end
    % NRMS = RMS / RMS noise
    for ss=1:length(sbj)
        for numCond=1:size(rmsSignal,2)
            nrms(ss,numCond) = rmsSignal(ss,numCond) / rmsNoise(ss,numCond) ;
        end
    end    
    % NRMSE = RMSE / RMS noise depending on the condition
    for ss=1:length(sbj)
        for numCond=2:4
            nrmse(ss,numCond-1) = rmse(ss,numCond-1) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,numCond)^2)); % divide by noise AM + noise from the condition
            nrmse(ss,numCond+3) = rmse(ss,numCond+3) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,numCond+5)^2)); % divide by noise AM + noise from the condition
        end
        % for s+t need to add both noise + AM
        nrmse(ss,4) = rmse(ss,4) / (sqrt(rmsNoise(ss,1)^2 + rmsNoise(ss,3)^2 + rmsNoise(ss,4)^2)); 
        nrmse(ss,8) = rmse(ss,8) / (sqrt(rmsNoise(ss,6)^2 + rmsNoise(ss,8)^2 + rmsNoise(ss,9)^2)); 
    end
    % NRMSE after regression for LR and SR (done only for aS+bT)
    load(['rmseCoefE' num2str(ee) '.mat']) % rmseS rmseT rmseS+T x LR/SR 
    load(['regCoefE' num2str(ee)]); % get the regression coef to get the amount of noise included
    for ss=1:length(sbj)
        nrmseCoef(ss,1) = rmseCoef(ss,3) / (sqrt( (coefF(ss,1,1)*rmsNoise(ss,1))^2 + (coefF(ss,1,2)*rmsNoise(ss,numCond))^2));
        nrmseCoef(ss,2) = rmseCoef(ss,6) / (sqrt( (coefF(ss,2,1)*rmsNoise(ss,6))^2 + (coefF(ss,2,2)*rmsNoise(ss,numCond+5))^2 )); 
    end
    figure('Renderer', 'painters', 'Position', [10 10 900 400])
    hold on;
    subplot(1,2,1)
    boxplot(1./[nrmse(:,1:4) nrmseCoef(:,1)])
    line([1 5],[1 1],'Color','r','LineWidth',2)
    title('LR')
    ylim([0 1.2])
    xticklabels({'lin','spat','temp','s+t','regress'})
    ylabel('1/NRMSE')
    subplot(1,2,2)
    boxplot(1./[nrmse(:,5:8) nrmseCoef(:,2)])
    title('SR')
    line([1 5],[1 1],'Color','r','LineWidth',2)
    xticklabels({'lin','spat','temp','s+t','regress'})
    ylabel('1/NRMSE')
    ylim([0 1.2])
    %%%%%% PB WITH THE REGRESSION BOXPLOT!!!!
    
    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Per electrode
    clear rmsNoiseChan rmsSignalChan rmseChan
    % RMS 
    for ss=1:length(sbj)
        for numCond=1:10
            for chan = 1 : size(sbj(ss,numCond).data.filteredWave,1)
                rmsNoiseChan(ss,numCond,chan) = rms( sbj(ss,numCond).data.noiseWave(chan,:) );
                rmsSignalChan(ss,numCond,chan) = rms( sbj(ss,numCond).data.filteredWave(chan,:) );
            end
        end
    end
    % RMSE 
    % spa, temp, s+t
    for ss=1:length(sbj)
        for numCond=2:5 
            rmseChan(ss,numCond-1,chan) = rms(sbj(ss,1).data.filteredWave(chan,:) -  sbj(ss,numCond).data.filteredWave(chan,:)); % long-range
            rmseChan(ss,numCond+3,chan) = rms(sbj(ss,6).data.filteredWave(chan,:) -  sbj(ss,numCond+5).data.filteredWave(chan,:)); % short-range
        end
    end    
    
end

