% regressions per subj
% same coef for all elec and time points
% linear set to 1

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
dataPath = {'/Users/Marlene/Documents/git/dataLR/LRlongDC/', '/Users/Marlene/Documents/git/dataLR/LRshortDC/'};


for ee=1:2
clear coefS coefT coefF
clear rmseSp rmseTe rmseFu rmseFuTopo correctedSpat correctedTemp subPredFul 
load([dataPath{ee} 'sbjprediction.mat'])

for ss=1:length(sbj)
    clear spatInt tempInt subSig predSp predTe predFul
    clear regCoefS regCoefT regCoefF  coefSp coefTe coefFu statsSp statsTe statsFu bintSp bintTe bintFu rSp rTe rFu rintSp rintTe rintFu
    
    
    %% get the spatial and temporal interaction terms
    % so substract the linear component from the prediction
    % also get the AM-linear that will be used in the regression
    for numCond=1:2 % long/short
        spatInt(:,:,numCond) = sbj(ss,3+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
        tempInt(:,:,numCond) = sbj(ss,4+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
        subSig(:,:,numCond) = sbj(ss,1+5*(numCond-1)).data.filteredWave - sbj(ss,2+5*(numCond-1)).data.filteredWave ;
    end
    
    % we want to set linear to 1 and only vary the coef for spat and
    % temp components
    % these regressions are therefore done on the residuals of the
    % difference AM-linear (the r2 are for the residuals as well)
    newDim = size(subSig,1)*size(subSig,2);
    for numCond=1:2 % long/short
        regCoefS = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1])]; % as 1*linear + spatial
        regCoefT = [ones(newDim,1) reshape(tempInt(:,:,numCond),[newDim 1])];
        regCoefF = [ones(newDim,1) reshape(spatInt(:,:,numCond),[newDim 1]) reshape(tempInt(:,:,numCond),[newDim 1])];
        [coefSp(numCond,:), bintSp(numCond,:,:), rSp(numCond,:), rintSp(numCond,:,:), statsSp(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefS);
        [coefTe(numCond,:), bintTe(numCond,:,:), rTe(numCond,:), rintTe(numCond,:,:), statsTe(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefT);
        [coefFu(numCond,:), bintFu(numCond,:,:), rFu(numCond,:), rintFu(numCond,:,:), statsFu(numCond,:)] = regress(reshape(subSig(:,:,numCond),[newDim 1]),regCoefF);
    end
    % vector STATS containing, in the following order, the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
    
    % calculate rmse after applying coef on the signal
    for numCond=1:2 % long/short
        predSp = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefSp(numCond,2)*spatInt(:,:,numCond);
        rmseSp(ss,numCond) = rms(predSp(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:) );
        
        predTe = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefTe(numCond,2)*tempInt(:,:,numCond);
        rmseTe(ss,numCond) = rms(predTe(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:) );
        
        predFul = sbj(ss,2+5*(numCond-1)).data.filteredWave + coefFu(numCond,2)*spatInt(:,:,numCond) + coefFu(numCond,3)*tempInt(:,:,numCond);
        rmseFu(ss,numCond) = rms( predFul(:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(:));
        % rms predFul for topo (per chan)
        for chan = 1: size(predFul,1)
            rmseFuTopo(ss,numCond,chan) = rms( predFul(chan,:) - sbj(ss,1+5*(numCond-1)).data.filteredWave(chan,:));
        end
     
        % new interaction terms with the regression terms
        correctedSpat(ss,numCond,:,:) = coefFu(numCond,2)*spatInt(:,:,numCond);
        correctedTemp(ss,numCond,:,:) = coefFu(numCond,3)*tempInt(:,:,numCond);
        subPredFul(ss,numCond,:,:) = predFul;
    end
    % regression coefficient for each sbj
    coefF(ss,:,:) = coefFu(:,2:3);
    coefS(ss,:,:) = coefSp(:,2);
    coefT(ss,:,:) = coefTe(:,2);
end

% save regression coefficients
save(['regCoefE' num2str(ee) '.mat'],'coefF','coefS','coefT');
% save rmse values for full regression on all electrodes (LR and SR)
save(['rmseFuTopo' num2str(ee) '.mat'],'rmseFuTopo');
% save rmse for spatial, temporal, s+t conditions x LR/SR
rmseCoef = [rmseSp(:,1) rmseTe(:,1) rmseFu(:,1) rmseSp(:,2) rmseTe(:,2) rmseFu(:,2)];
save(['rmseCoefE' num2str(ee) '.mat'],'rmseCoef'); 
% save spatial and temporal regression prediction
save(['subPredFulE' num2str(ee) '.mat'],'subPredFul'); 

end
