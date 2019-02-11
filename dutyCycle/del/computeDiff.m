function [difference] = computeDiff(data1,data2)
% input 2 Axx
% output the difference

%%% Difference
% to compute the difference, can do just the substraction for the wave, 
% but for the amplitude, needs to be computed from the sin and cos
% difference

difference.wave = data1.wave-data2.wave;

difference.sin = data1.sin - data2.sin;
difference.cos = data1.cos - data2.cos;
difference.amp = sqrt(difference.sin.^2 + difference.cos.^2);


% values that do not change
difference.nt = data1.nt;
difference.nfr = data1.nfr;
difference.time = data1.time;
difference.pval = NaN(size(data1.pval));
difference.confradius = NaN(size(data1.confradius));
difference.i1f1 = data1.i1f1;
difference.elec = data1.elec;
difference.dtms = data1.dtms;
difference.nchan = data1.nchan;
difference.ndft = data1.ndft;
difference.freq = data1.freq;
difference.label = data1.label;

end