%%%
% could add the signal from the 2 s epochs (cleanData) or from the 400 ms
% cycle (Axx). To keep the structure of the data, do it from the Axx

% 2 ways to do the linear summation
% sum the waves and then do fft
% the pb is that because the epoch is 400 ms, will loose many frequencies
% do the summation of the cos and sin
% but the shift for computing the motion is tricky to do

load AxxS01_V2.mat
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
elect = 39;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULTANEOUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 is the original long range
% 2 the prediction long range
% 3 is the original short range
% 4 is the prediction short range
% 5 the difference long range (original - prediction) -> take the absolute
% value for the amp otherwise won't plot the negative 
% 6 the difference short range

Simult(1) = Axx(4);
Simult(3) = Axx(8);

% stuff needed but does not change
Simult(2).time =  Axx(1).time ;
Simult(2).nt =  Axx(1).nt;
Simult(2).freq =  Axx(1).freq; % in this case, same frequencies
% have to find out later, for now just put NaN
Simult(2).pval =  NaN(size(Axx(1).pval));
Simult(2).confradius =  NaN(size(Axx(1).confradius));
% do the addition
Simult(2).wave = Axx(2).wave + Axx(3).wave;
Simult(2).sin = Axx(2).sin + Axx(3).sin;
Simult(2).cos = Axx(2).cos + Axx(3).cos;
Simult(2).amp = sqrt(Simult(2).sin.^2 + Simult(2).cos.^2);


% do the addition
Simult(4).wave = Axx(6).wave + Axx(7).wave;
Simult(4).sin = Axx(6).sin + Axx(7).sin;
Simult(4).cos = Axx(6).cos + Axx(7).cos;
Simult(4).amp = sqrt(Simult(4).sin.^2 + Simult(4).cos.^2);


% stuff needed for plotting but does not change
for nb=1:6
    Simult(nb).time =  Axx(1).time ;
    Simult(nb).nt =  Axx(1).nt;
    Simult(nb).freq =  Axx(1).freq; % in this case, same frequencies
    Simult(nb).pval =  NaN(size(Axx(1).pval));
    Simult(nb).confradius =  NaN(size(Axx(1).confradius));
end

Simult(5).wave = Simult(1).wave-Simult(2).wave;
Simult(6).wave = Simult(3).wave-Simult(4).wave;
Simult(5).amp = abs(Simult(1).amp-Simult(2).amp);
Simult(6).amp = abs(Simult(3).amp-Simult(4).amp);
Simult(5).sin = Simult(1).sin-Simult(2).sin;
Simult(6).sin = Simult(3).sin-Simult(4).sin;
Simult(5).cos = Simult(1).cos-Simult(2).cos;
Simult(6).cos = Simult(3).cos-Simult(4).cos;

interactiveSteadyStatePlot_noY(cfg,Simult)


% plot difference for the short and long range
% this is just picking one electrode: 
elect = 39; % PO8
diffLR = Simult(1).wave-Simult(2).wave;
diffSR = Simult(3).wave-Simult(4).wave;
diffLRamp = Simult(1).amp-Simult(2).amp;
diffSRamp = Simult(3).amp-Simult(4).amp;
figure; subplot(2,1,1);plot(Simult(1).time,diffLR(elect,:)); hold on; plot(Simult(1).time,diffSR(elect,:));legend('LR','SR')
subplot(2,1,2);
plot(Simult(1).freq(Simult(1).freq<30),diffLRamp(elect,Simult(1).freq<30)); hold on;
plot(Simult(1).freq(Simult(1).freq<30),diffSRamp(elect,Simult(1).freq<30));

% plot(Simult(4).wave(elect,:)); hold on;
% plot(Axx(6).wave(elect,:),'r'); 
% plot(Axx(7).wave(elect,:),'g'); 
% plot(Axx(8).wave(elect,:),'k'); 











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% attention in this case Axx(3) has to be shifted in time
Motion(1) = Axx(1);

% compute the new Axx(3) with a shift in time of 200 ms (half cycle)
% first for the wave
shiftLR.wave = circshift(Axx(3).wave,[0 length(Axx(3).wave)/2]);
% then for the cos and sin
% for this needs the phase shift calculated for each frequency
% phase shift = shift (200 ms) / periods * 360 degrees
periods = 1 ./ Axx(3).freq;
phaseShift = 0.200 ./ periods .* 360;

% trigonometry.. 
% cos' = A*cos(theta + phi) = cos(theta)*alpha*cos(phi)
% (theta = angle difference, phi = original angle)
% cos' = cosd(phaseShift) * orignialCos - sind(phaseShift) * orignialSin
% Similar stuff for sin' but different "identification"
% sin' = sind(phaseShift) * originalCos + cosd(phaseShift) * originalSin
shiftLR.sin = sind(phaseShift) .* Axx(3).cos + cosd(phaseShift) .* Axx(3).sin;
shiftLR.cos = cosd(phaseShift) .* Axx(3).cos - sind(phaseShift) .* Axx(3).sin;

% now do the addition
Motion(2).wave = Axx(2).wave +shiftLR.wave;
Motion(2).sin = Axx(2).sin + shiftLR.sin;
Motion(2).cos = Axx(2).cos + shiftLR.cos;
Motion(2).amp = sqrt(Motion(2).sin.^2 + Motion(2).cos.^2);
Motion(2).freq =  Axx(1).freq; % in this case, same frequencies


% do the same for the SR
Motion(3) = Axx(5);

shiftSR.wave = circshift(Axx(7).wave,[0 length(Axx(7).wave)/2]);
periods = 1 ./ Axx(7).freq;
phaseShift = 0.200 ./ periods .* 360;
shiftSR.sin = sind(phaseShift) .* Axx(7).cos + cosd(phaseShift) .* Axx(7).sin;
shiftSR.cos = cosd(phaseShift) .* Axx(7).cos - sind(phaseShift) .* Axx(7).sin;

% now do the addition
Motion(4).wave = Axx(6).wave +shiftSR.wave;
Motion(4).sin = Axx(6).sin + shiftSR.sin;
Motion(4).cos = Axx(6).cos + shiftSR.cos;
Motion(4).amp = sqrt(Motion(4).sin.^2 + Motion(4).cos.^2);
Motion(4).freq =  Axx(1).freq; % in this case, same frequencies



% stuff needed for plotting but does not change
for nb=1:6
    Motion(nb).time =  Axx(1).time ;
    Motion(nb).nt =  Axx(1).nt;
    Motion(nb).freq =  Axx(1).freq; 
    Motion(nb).pval =  NaN(size(Axx(1).pval));
    Motion(nb).confradius =  NaN(size(Axx(1).confradius));
end

Motion(5).wave = Motion(1).wave-Motion(2).wave;
Motion(6).wave = Motion(3).wave-Motion(4).wave;
Motion(5).amp = abs(Motion(1).amp-Motion(2).amp);
Motion(6).amp = abs(Motion(3).amp-Motion(4).amp);
Motion(5).sin = Motion(1).sin-Motion(2).sin;
Motion(6).sin = Motion(3).sin-Motion(4).sin;
Motion(5).cos = Motion(1).cos-Motion(2).cos;
Motion(6).cos = Motion(3).cos-Motion(4).cos;

interactiveSteadyStatePlot(cfg,Motion)

% plot difference for the short and long range
motionDiffLR = Motion(1).wave-Motion(2).wave;
motionDiffSR = Motion(3).wave-Motion(4).wave;
motionDiffLRamp = Motion(1).amp-Motion(2).amp;
motionDiffSRamp = Motion(3).amp-Motion(4).amp;
figure; subplot(2,1,1);plot(Motion(1).time,motionDiffLR(elect,:)); hold on; plot(Motion(1).time,motionDiffSR(elect,:));legend('LR','SR')
subplot(2,1,2);
plot(Motion(1).freq(Motion(1).freq<30),motionDiffLRamp(elect,Motion(1).freq<30)); hold on;
plot(Motion(1).freq(Motion(1).freq<30),motionDiffSRamp(elect,Motion(1).freq<30));





% check that the results are the same when doing fft on the sum of the 2
% signals instead of addition of cos and sin
checkSum.wave = Axx(2).wave +shiftLR.wave;
ndft = length(shiftLR.wave);
nfr = floor(ndft/2)+1;
dft = dftmtx(ndft); % should be 204*204
dftData = checkSum.wave*dft; %Do the fourier transform of the data. should be 204
dftData = dftData(:,1:nfr); %Select just the unique frequencies. should be 103 
freqsToDouble = 2:(nfr-1+mod(ndft,2)); 
dftData(:,freqsToDouble) = 2*dftData(:,freqsToDouble);
dftData = dftData/ndft;   
for iChan = 1:Axx(2).nchan    
    checkSum.amp(iChan,:) = abs(dftData(iChan,:));
    checkSum.cos(iChan,:) = real(dftData(iChan,:));
    checkSum.sin(iChan,:) = -imag(dftData(iChan,:));
end
checkSum.sample  = Axx(1).fsample ;
checkSum.freq = (checkSum.sample/2)*linspace(0,1,nfr);
checkSum.time =  Axx(1).time ;
checkSum.nt =  Axx(1).nt;
checkSum.pval =  NaN(size(Axx(1).pval));
checkSum.confradius =  NaN(size(Axx(1).confradius));

interactiveSteadyStatePlot(cfg,checkSum)
interactiveSteadyStatePlot(cfg,Motion)

for dp=1:103
    test.amp(:,dp) = Motion(2).amp(:,(dp-1)*5+1);
end
figure;plot(Motion(2).freq(1:40),Motion(2).amp(39,1:40))
figure;plot(checkSum.freq(1:20),checkSum.amp(39,1:20))
