function [Motion] = predictMotion(dataAxx)
% 1 is the original long range
% 2 the prediction long range
% 3 is the original short range
% 4 is the prediction short range
% 5 the difference long range (original - prediction) -> take the absolute
% value for the amp otherwise won't plot the negative 
% 6 the difference short range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% attention dataAxx(3) has to be shifted in time
Motion(1) = dataAxx(1);

% compute the new dataAxx(3) with a shift in time of 200 ms (half cycle)
% first for the wave
shiftLR.wave = circshift(dataAxx(3).wave,[0 length(dataAxx(3).wave)/2]);
% then for the cos and sin
% for this needs the phase shift calculated for each frequency
% phase shift = shift (200 ms) / periods * 360 degrees
periods = 1 ./ dataAxx(3).freq;
phaseShift = 0.200 ./ periods .* 360;

% trigonometry.. 
% cos' = A*cos(theta + phi) = cos(theta)*alpha*cos(phi)
% (theta = angle difference, phi = original angle)
% cos' = cosd(phaseShift) * orignialCos - sind(phaseShift) * orignialSin
% Similar stuff for sin' but different "identification"
% sin' = sind(phaseShift) * originalCos + cosd(phaseShift) * originalSin
shiftLR.sin = sind(phaseShift) .* dataAxx(3).cos + cosd(phaseShift) .* dataAxx(3).sin;
shiftLR.cos = cosd(phaseShift) .* dataAxx(3).cos - sind(phaseShift) .* dataAxx(3).sin;

% now do the addition
Motion(2).wave = dataAxx(2).wave +shiftLR.wave;
Motion(2).sin = dataAxx(2).sin + shiftLR.sin;
Motion(2).cos = dataAxx(2).cos + shiftLR.cos;
Motion(2).amp = sqrt(Motion(2).sin.^2 + Motion(2).cos.^2);
Motion(2).freq =  dataAxx(1).freq; % in this case, same frequencies


% do the same for the SR
Motion(3) = dataAxx(7);

shiftSR.wave = circshift(dataAxx(9).wave,[0 length(dataAxx(9).wave)/2]);
periods = 1 ./ dataAxx(9).freq;
phaseShift = 0.200 ./ periods .* 360;
shiftSR.sin = sind(phaseShift) .* dataAxx(9).cos + cosd(phaseShift) .* dataAxx(9).sin;
shiftSR.cos = cosd(phaseShift) .* dataAxx(9).cos - sind(phaseShift) .* dataAxx(9).sin;

% now do the addition
Motion(4).wave = dataAxx(8).wave +shiftSR.wave;
Motion(4).sin = dataAxx(8).sin + shiftSR.sin;
Motion(4).cos = dataAxx(8).cos + shiftSR.cos;
Motion(4).amp = sqrt(Motion(4).sin.^2 + Motion(4).cos.^2);
Motion(4).freq =  dataAxx(1).freq; % in this case, same frequencies


%%% Difference
% to compute the difference, can do just the substraction for the wave, 
% but for the amplitude, needs to be computed from the sin and cos
% difference
Motion(5).wave = Motion(1).wave-Motion(2).wave;
Motion(6).wave = Motion(3).wave-Motion(4).wave;

Motion(5).sin = Motion(1).sin - Motion(2).sin;
Motion(6).sin = Motion(3).sin  - Motion(4).sin;
Motion(5).cos = Motion(1).cos - Motion(2).cos;
Motion(6).cos = Motion(3).cos - Motion(4).cos;
Motion(5).amp = sqrt(Motion(5).sin.^2 + Motion(5).cos.^2);
Motion(6).amp = sqrt(Motion(6).sin.^2 + Motion(6).cos.^2);


% stuff needed for plotting but does not change
for nb=1:6
    Motion(nb).time =  dataAxx(1).time ;
    Motion(nb).nt =  dataAxx(1).nt;
    Motion(nb).freq =  dataAxx(1).freq; 
    Motion(nb).pval =  NaN(size(dataAxx(1).pval));
    Motion(nb).confradius =  NaN(size(dataAxx(1).confradius));
end


%%% add labels
% 1 is the original long range
Motion(1).label = 'original LR';
% 2 the prediction long range
Motion(2).label = 'prediction LR';
% 3 is the original short range
Motion(3).label = 'original SR';
% 4 is the prediction short range
Motion(4).label = 'prediction SR';
% 5 the difference long range (original - prediction) 
Motion(5).label = 'difference LR';
% 6 the difference short range
Motion(6).label = 'difference SR';

end

% % check that the results are the same when doing fft on the sum of the 2
% % signals instead of addition of cos and sin
% checkSum.wave = dataAxx(2).wave +shiftLR.wave;
% ndft = length(shiftLR.wave);
% nfr = floor(ndft/2)+1;
% dft = dftmtx(ndft); % should be 204*204
% dftData = checkSum.wave*dft; %Do the fourier transform of the data. should be 204
% dftData = dftData(:,1:nfr); %Select just the unique frequencies. should be 103 
% freqsToDouble = 2:(nfr-1+mod(ndft,2)); 
% dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
% dftData = dftData/ndft;   
% for iChan = 1:dataAxx(2).nchan    
%     checkSum.amp(iChan,:) = abs(dftData(iChan,:));
%     checkSum.cos(iChan,:) = real(dftData(iChan,:));
%     checkSum.sin(iChan,:) = -imag(dftData(iChan,:));
% end
% checkSum.sample  = dataAxx(1).fsample ;
% checkSum.freq = (checkSum.sample/2)*linspace(0,1,nfr);
% checkSum.time =  dataAxx(1).time ;
% checkSum.nt =  dataAxx(1).nt;
% checkSum.pval =  NaN(size(dataAxx(1).pval));
% checkSum.confradius =  NaN(size(dataAxx(1).confradius));
% 
% interactiveSteadyStatePlot(cfg,checkSum)
% interactiveSteadyStatePlot(cfg,Motion)
% 
% for dp=1:103
%     test.amp(:,dp) = Motion(2).amp(:,(dp-1)*5+1);
% end
% figure;plot(Motion(2).freq(1:40),Motion(2).amp(39,1:40))
% figure;plot(checkSum.freq(1:20),checkSum.amp(39,1:20))
