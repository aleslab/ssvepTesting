close all;
clearvars;
% to get an approximation of the noise level, need to compare the signal
% 
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults

dataPath = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};

for dd=1:2 % which experiment (different duty-cycle used)
    load([dataPath{dd} 'sbjprediction.mat'])
    
    for sbjInd=1:length(sbj)
        % only keep neighbouring freq
        for condIdx=1:size(sbj,2)
            filtIdx = determineFilterIndices( 'nf1low50', sbj(sbjInd,condIdx).data.freq, sbj(sbjInd,condIdx).data.i1f1 );
            filtVoisin = [filtIdx-1 filtIdx+1 ];
            
            %Create a logical matrix selecting frequency components.
            filtMat = false(size(sbj(sbjInd,condIdx).data.amp));
            filtMat(:,filtVoisin) = true;
            
            %Combine the filter and sig vaules with a logical AND.
            % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
            sbj(sbjInd,condIdx).data.activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
            cfg.activeFreq =  sbj(sbjInd,condIdx).data.activeFreq;
            
            [ sbj(sbjInd,condIdx).data.baselineWave ] = filterSteadyState( cfg, sbj(sbjInd,condIdx).data );
        end
        % all freq but harmonics
        for condIdx=1:size(sbj,2)
            filtIdx = determineFilterIndices( 'nf1low50', sbj(sbjInd,condIdx).data.freq, sbj(sbjInd,condIdx).data.i1f1 );
            
            %Create a logical matrix selecting frequency components.
            filtMat = true(size(sbj(sbjInd,condIdx).data.amp));
            filtMat(:,filtIdx) = false;
            filtMat(:,filtIdx(end):end) = false;
            
            %Combine the filter and sig vaules with a logical AND.
            % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
            sbj(sbjInd,condIdx).data.activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
            cfg.activeFreq =  sbj(sbjInd,condIdx).data.activeFreq;
            
            [ sbj(sbjInd,condIdx).data.baselineWave ] = filterSteadyState( cfg, sbj(sbjInd,condIdx).data );
        end
        %%% Correction shift cycle only for E1
        if dd==1
            for condIx = 1:length(sbj)
                sbj(sbjInd,condIdx).data.baselineWave = circshift(sbj(sbjInd,condIdx).data.baselineWave,[0 length(sbj(sbjInd,condIdx).data.baselineWave)/2]);
            end
        end
    end
    
    save([dataPath{dd} 'sbjprediction.mat'],'sbj');
end





