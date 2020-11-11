
function Axx = reformatAxx(AxxOld)
% reformat Axx so that it works with other programs

oldNames = fieldnames(AxxOld);
newNames = lower(oldNames);
for k=1:length(oldNames)
    Axx.(newNames{k}) = AxxOld.(oldNames{k}) ;
end
Axx.time = linspace(0,(Axx.nt-1)*Axx.dtms,Axx.nt);
Axx.pval = zeros(128,108);
Axx.amp = Axx.amp';
Axx.cos = Axx.cos';
Axx.sin = Axx.sin';
Axx.wave = Axx.wave';
Axx.confradius = Axx.pval;
Axx.freq = linspace(0,(Axx.nfr-1)*Axx.dfhz,Axx.nfr)';

end