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
saveplot = 0;

pred = []; intera = [];

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
    
    pred = [pred avPredictions(1:10)];
    intera = [intera avInteractions(1:8)];
end

save('results.mat','pred')









%%%% to get SSVEP cannot just do a substraction (as with wave)
%%%%%% compute differences
pickCond = [2:5 ; 7:10; 12:15; 17:20];
for cc=1:size(pickCond,1)
    condCompare = pickCond(cc,:);
    fixCond = condCompare(1) - 1;
    for cond=1:length(condCompare)
        diff(cc,cond) = computeDiff(pred(fixCond),pred(condCompare(cond)));
    end
end
% create a new field "filteredWave" after low-pass filter
% warning about filter is NORMAL
for cc=1:4
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

save('diff.mat','diff')
