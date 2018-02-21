function [Motion] = reconstructTiming(AxxToShift,AxxNoChange)
% input 2 Axx files, 1st is the one to shift in time, the 2nd is not
% shifted, then the 2 signals are summed
% output the sum of the 2 signals

% compute the new Axx with a shift in time of 200 ms (half cycle)
% first for the wave
shift.wave = circshift(AxxToShift.wave,[0 length(AxxToShift.wave)/2]);
% then for the cos and sin
% for this needs the phase shift calculated for each frequency
% phase shift = shift (200 ms) / periods * 360 degrees
periods = 1 ./ AxxToShift.freq;
phaseShift = 0.2 ./ periods .* 360;

% trigonometry.. 
% cos' = A*cos(theta + phi) = cos(theta)*alpha*cos(phi)
% (theta = angle difference, phi = original angle)
% cos' = cosd(phaseShift) * orignialCos - sind(phaseShift) * orignialSin
% Similar stuff for sin' but different "identification"
% sin' = sind(phaseShift) * originalCos + cosd(phaseShift) * originalSin
shift.sin = sind(phaseShift) .* AxxToShift.cos + cosd(phaseShift) .* AxxToShift.sin;
shift.cos = cosd(phaseShift) .* AxxToShift.cos - sind(phaseShift) .* AxxToShift.sin;

% now do the addition
Motion.wave = AxxNoChange.wave +shift.wave;
Motion.sin = AxxNoChange.sin + shift.sin;
Motion.cos = AxxNoChange.cos + shift.cos;
Motion.amp = sqrt(Motion.sin.^2 + Motion.cos.^2);
Motion.freq =  AxxNoChange.freq; % in this case, same frequencies


% values that do not change
Motion.nt = AxxNoChange.nt;
Motion.nfr = AxxNoChange.nfr;
Motion.time = AxxNoChange.time;
Motion.pval = NaN(size(AxxNoChange.pval));
Motion.confradius = NaN(size(AxxNoChange.confradius));
Motion.i1f1 = AxxNoChange.i1f1;
Motion.elec = AxxNoChange.elec;
Motion.dtms = AxxNoChange.dtms;
Motion.nchan = AxxNoChange.nchan;
Motion.ndft = AxxNoChange.ndft;
Motion.freq = AxxNoChange.freq;
Motion.label = AxxNoChange.label;

end