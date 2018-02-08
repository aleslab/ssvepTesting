function [groupAvMotion, groupAvSimult] = doAverage(dataDir,listData,cond)

for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    dataWave(:,:,ff) = Axx.predictMotion(cond).wave;
    dataAmp(:,:,ff) = Axx.predictMotion(cond).amp;
    dataSin(:,:,ff) = Axx.predictMotion(cond).sin;
    dataCos(:,:,ff) = Axx.predictMotion(cond).cos;
    dataWaveSimult(:,:,ff) = Axx.predictSimult(cond).wave;
    dataAmpSimult(:,:,ff) = Axx.predictSimult(cond).amp;
    dataSinSimult(:,:,ff) = Axx.predictSimult(cond).sin;
    dataCosSimult(:,:,ff) = Axx.predictSimult(cond).cos;
end

% values that do not change
groupAvMotion.nt = Axx.predictMotion(1).nt;
groupAvMotion.nfr = Axx.predictMotion(1).nfr;
groupAvMotion.time = Axx.predictMotion(1).time;
groupAvMotion.pval = NaN(size(Axx.predictMotion(1).pval));
groupAvMotion.confradius = NaN(size(Axx.predictMotion(1).confradius));
groupAvMotion.i1f1 = Axx.predictMotion(1).i1f1;
groupAvMotion.elec = Axx.predictMotion(1).elec;
groupAvMotion.dtms = Axx.predictMotion(1).dtms;
groupAvMotion.nchan = Axx.predictMotion(1).nchan;
groupAvMotion.ndft = Axx.predictMotion(1).ndft;

% copy in the other average
groupAvSimult = groupAvMotion;

% depends on the condition so do it separately
groupAvMotion.wave = mean(dataWave,3);
groupAvMotion.amp = mean(dataAmp,3);
groupAvMotion.sin = mean(dataSin,3);
groupAvMotion.cos = mean(dataCos,3);
groupAvMotion.freq = Axx.predictMotion(cond).freq;
groupAvMotion.label = Axx.predictMotion(cond).label;

groupAvSimult.wave = mean(dataWaveSimult,3);
groupAvSimult.amp = mean(dataAmpSimult,3);
groupAvSimult.sin = mean(dataSinSimult,3);
groupAvSimult.cos = mean(dataCosSimult,3);
groupAvSimult.freq = Axx.predictSimult(cond).freq;
groupAvSimult.label = Axx.predictSimult(cond).label;

end




