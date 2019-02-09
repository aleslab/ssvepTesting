close all;
clearvars;
% create a filteredWave Field

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};

for dd=1:2 % which experiment (different duty-cycle used)
    load([dataPath{dd} 'sbjprediction.mat'])

    for sbjInd=1:length(sbj)
    for condIdx=1:size(sbj,2)
        filtIdx = determineFilterIndices( 'nf1low50', sbj(sbjInd,condIdx).data.freq, sbj(sbjInd,condIdx).data.i1f1 );
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(sbj(sbjInd,condIdx).data.amp));
        filtMat(:,filtIdx) = true;
        
        %Combine the filter and sig vaules with a logical AND.
        % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        sbj(sbjInd,condIdx).data.activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
        cfg.activeFreq =  sbj(sbjInd,condIdx).data.activeFreq;
        
        [ sbj(sbjInd,condIdx).data.filteredWave ] = filterSteadyState( cfg, sbj(sbjInd,condIdx).data );
    end    
    
    %%% Correction shift cycle only for E1
    if dd==1
        for condIx = 1:length(sbj)
            sbj(sbjInd,condIdx).data.filteredWave = circshift(sbj(sbjInd,condIdx).data.filteredWave,[0 length(sbj(sbjInd,condIdx).data.filteredWave)/2]);
        end
    end
    end
    
    save([dataPath{dd} 'sbjprediction.mat'],'sbj');
end


