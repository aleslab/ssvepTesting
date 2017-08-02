function [barH sigH] = pdSpecPlot(freqs,amps,sig,varargin)
%function pdSpecPlot(freqs,amps,sig,varargin)


nFreqs = length(freqs);

barcolor = 'r';
plotErrorBar = false;

%Makre sure data is a 1xn vector. 
amps = squeeze(amps);
dataSz = size(amps);
if dataSz(2)<dataSz(1)
    amps = amps';
end


if ~isempty(varargin)
    
    for iArg = 1:2:length(varargin)

        thisParam = lower(varargin{iArg});

        switch thisParam

            case 'errorbar'
                plotErrorBar = varargin{iArg+1};
                
            case 'snr'
                
                if varargin{iArg+1},

                    noiseList = setdiff(2:nFreqs,sig);
                    
                    idx = 1;
                    for iNoise = noiseList,
                        
                        if ismember(iNoise,sig) || ismember(iNoise+1,sig)
                            continue;
                        end
                        
                        if iNoise==1  
                            noiseLevelMean(idx) = mean(amps(iNoise:iNoise+1));
                        elseif iNoise == nFreqs
                            noiseLevelMean(idx) = mean(amps(iNoise-1:iNoise));
                        else
                            noiseLevelMean(idx) = mean( amps([iNoise iNoise+1]) );
                        end
                        idx = idx+1;
                    end
                    
                            
                    resampIdx = round(linspace(1,length(noiseLevelMean),nFreqs));
                    
                    amps = (amps./noiseLevelMean(resampIdx));
%                   amps = amps-noiseLevelMean(resampIdx);
                    
                end

            case 'color'
                barcolor = varargin{iArg+1};
                                        
            otherwise 
                disp(['Unknown paramter ' thisParam])
            
        end
        
    end
else
    plotErrorBar = false;
end

barWidth = .8;
barH = bar(freqs,amps,barWidth,'w');
hold on
set(barH,'edgecolor','k');

if max(amps)~=0 && min(freqs)<40
    
axis([freqs(1) 40 0 max(amps(2:end))*1.1 ]);

end

colorVec = zeros(size(amps));

colorVec(sig) = amps(sig);


sigH = bar(freqs,colorVec,barWidth);


if plotErrorBar  
    %Noise estimate
    noise1 = amps(sig+1);
    noise2 = amps(sig-1);
    noiseEst = (noise1+noise2)/2;
    errorbar(freqs(sig),amps(sig),zeros(size(noiseEst)),noiseEst,'.k');
end

set(sigH,'facecolor',barcolor);
hold off

ylabel('uV')
xlabel('Frequency (Hz)')


