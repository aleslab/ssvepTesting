close all;
clearvars;

% to get an approximation of baseline brain activity, get the amount of
% noise in the harmonics compare to sig
% to approximate this baseline noise/activity, get the  average
% amplitude of the neighbouring freq harmonic


            
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
            for chan=1:size(sbj(sbjInd,condIdx).data.amp,1)
                baseAmp(chan,:) = mean([sbj(sbjInd,condIdx).data.amp(chan,filtIdx-1); sbj(sbjInd,condIdx).data.amp(chan,filtIdx+1)]);
            end
            
    % reconstruct waveform
    % this is actually a WRONG waveform as we cannot know the phase
    % of the signal. However, it can be used to calculate RMS or
    % sumofsquares because the deviation over time will be the same
    % whatever the phase
    
    %3 rebuild original source signal
amp = abs(ya_fft);
phase = unwrap(angle(ya_fft));
ya_newifft=ifft(mag.*exp(i*phase));
    steadystate.amp(iChan,:) = abs(meanDftData);
    steadystate.cos(iChan,:) = real(meanDftData);
    steadystate.sin(iChan,:) = -imag(meanDftData);
    %% How do I reconstruct waveform with only the amplitudes?? 
    % inverse fourier transform 
    % keep the cos and sin of the original signal (although it is not
    % correct)
    waveRecon = real(dft(:,selFr))*sbj(sbjInd,condIdx).data.cos(iChan,selFr)' - imag(dft(:,selFr))*sbj(sbjInd,condIdx).data.sin(iChan,selFr)';    

        end
    end
    
%     for sbjInd=1:length(sbj)
%         % only keep neighbouring freq
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
%         % all freq but harmonics
%         for condIdx=1:size(sbj,2)
%             filtIdx = determineFilterIndices( 'nf1low50', sbj(sbjInd,condIdx).data.freq, sbj(sbjInd,condIdx).data.i1f1 );
%             
%             %Create a logical matrix selecting frequency components.
%             filtMat = true(size(sbj(sbjInd,condIdx).data.amp));
%             filtMat(:,filtIdx) = false;
%             filtMat(:,filtIdx(end):end) = false;
%             
%             %Combine the filter and sig vaules with a logical AND.
%             % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
%             sbj(sbjInd,condIdx).data.activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
%             cfg.activeFreq =  sbj(sbjInd,condIdx).data.activeFreq;
%             
%             [ sbj(sbjInd,condIdx).data.baselineWave ] = filterSteadyState( cfg, sbj(sbjInd,condIdx).data );
%         end
%         %%% Correction shift cycle only for E1
%         if dd==1
%             for condIx = 1:length(sbj)
%                 sbj(sbjInd,condIdx).data.baselineWave = circshift(sbj(sbjInd,condIdx).data.baselineWave,[0 length(sbj(sbjInd,condIdx).data.baselineWave)/2]);
%             end
%         end
%     end
%     
%     save([dataPath{dd} 'sbjprediction.mat'],'sbj');
end

dft = dftmtx(sbj(1,1).data.ndft);
dft = dft(:,1:sbj(1,1).data.nfr); 

sbj(1,1).data.amp(23,7)
sbj(1,1).data.cos(23,7)
sbj(1,1).data.sin(23,7)
test=real(dft(:,7))*sbj(1,1).data.cos(23,7)/2' - imag(dft(:,7))*sbj(1,1).data.sin(23,7)/2'  ;  
test=reshape(test,sbj(1,1).data.nt,[])';
figure;plot(mean(test))
    
sbj(1,1).data.amp(23,6)
sbj(1,1).data.cos(23,6)
sbj(1,1).data.sin(23,6)
testN=real(dft(:,6))*sbj(1,1).data.cos(23,6)' - imag(dft(:,6))*sbj(1,1).data.sin(23,6)'  ;  
testN=reshape(testN,sbj(1,1).data.nt,[])';
hold on; 
plot(mean(testN))

