function [groupAv] = doAverage(dataDir,listData,cond)

for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    dataWave(:,:,ff) = Axx.predictMotion(cond).wave;
    dataAmp(:,:,ff) = Axx.predictMotion(cond).amp;
    dataSin(:,:,ff) = Axx.predictMotion(cond).sin;
    dataCos(:,:,ff) = Axx.predictMotion(cond).cos;
end

groupAv.wave = mean(dataWave,3);
groupAv.amp = mean(dataAmp,3);
groupAv.sin = mean(dataSin,3);
groupAv.cos = mean(dataCos,3);
groupAv.freq = Axx.predictMotion(cond).freq;
groupAv.nt = Axx.predictMotion(1).nt;
groupAv.nfr = Axx.predictMotion(1).nfr;
groupAv.time = Axx.predictMotion(1).time;
groupAv.pval = NaN(size(Axx.predictMotion(1).pval));
groupAv.confradius = NaN(size(Axx.predictMotion(1).confradius));
groupAv.label = Axx.predictMotion(cond).label;
groupAv.i1f1 = Axx.predictMotion(1).i1f1;
groupAv.elec = Axx.predictMotion(1).elec;
groupAv.dtms = Axx.predictMotion(1).dtms;
groupAv.nchan = Axx.predictMotion(1).nchan;
groupAv.ndft = Axx.predictMotion(1).ndft;
end


