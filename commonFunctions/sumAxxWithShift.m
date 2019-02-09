function [motion] = sumAxxWithShift(AxxToShift,AxxNoChange,reshaping)
% input 2 Axx files, 1st is the one to shift in time, the 2nd is not
% shifted, then the 2 signals are summed
% output the sum of the 2 signals
% reshaping = optional argument = new reshape time window if the output needs to be reshaped to a new time
% window (cycle)

% compute the new Axx with a shift in time of 200 ms (half cycle)
% first for the wave
shift.wave = circshift(AxxToShift.wave,[0 length(AxxToShift.wave)/2]);
% then for the cos and sin
% for this needs the phase shift calculated for each frequency
% phase shift = shift (200 ms) / periods * 360 degrees
periods = 1 ./ AxxToShift.freq;
phaseShift = AxxToShift.time(end/2+1)/1000 ./ periods .* 360; 

% trigonometry.. 
% cos' = A*cos(theta + phi) = cos(theta)*alpha*cos(phi)
% (theta = angle difference, phi = original angle)
% cos' = cosd(phaseShift) * orignialCos - sind(phaseShift) * orignialSin
% Similar stuff for sin' but different "identification"
% sin' = sind(phaseShift) * originalCos + cosd(phaseShift) * originalSin
shift.sin = sind(phaseShift) .* AxxToShift.cos + cosd(phaseShift) .* AxxToShift.sin;
shift.cos = cosd(phaseShift) .* AxxToShift.cos - sind(phaseShift) .* AxxToShift.sin;

% now do the addition
motion.wave = AxxNoChange.wave +shift.wave;
motion.sin = AxxNoChange.sin + shift.sin;
motion.cos = AxxNoChange.cos + shift.cos;
motion.amp = sqrt(motion.sin.^2 + motion.cos.^2);

if nargin == 3 
    for ch=1:size(motion.wave,1)
        shortWave(ch,:,:) = reshape(motion.wave(ch,:),reshaping,[])';
    end
    motion.wave = squeeze(mean(shortWave,2));
    motion.time = AxxNoChange.time(1:reshaping);
else
    motion.time = AxxNoChange.time;
end

% values that do not change
motion.nt = AxxNoChange.nt;
motion.nfr = AxxNoChange.nfr;
motion.time = AxxNoChange.time;
motion.pval = NaN(size(AxxNoChange.pval));
motion.confradius = NaN(size(AxxNoChange.confradius));
motion.i1f1 = AxxNoChange.i1f1;
motion.elec = AxxNoChange.elec;
motion.dtms = AxxNoChange.dtms;
motion.nchan = AxxNoChange.nchan;
motion.ndft = AxxNoChange.ndft;
motion.freq = AxxNoChange.freq;
motion.label = AxxNoChange.label;

end