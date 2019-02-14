function [groupAv] = averageSbj(dataSbj)

for ff=1:size(dataSbj,1)
    for cond=1:size(dataSbj,2)
        dataWave(:,:,cond,ff) = dataSbj(ff,cond).data.wave;
        dataAmp(:,:,cond,ff) = dataSbj(ff,cond).data.amp;
        dataSin(:,:,cond,ff) = dataSbj(ff,cond).data.sin;
        dataCos(:,:,cond,ff) = dataSbj(ff,cond).data.cos;
        dataFil(:,:,cond,ff) = dataSbj(ff,cond).data.filteredWave;
        dataNoise(:,:,cond,ff) = dataSbj(ff,cond).data.noiseWave;
    end
end

for cond=1:size(dataSbj,2) 
% values that do not change
groupAv(cond).nt = dataSbj(ff,cond).data.nt;
groupAv(cond).nfr = dataSbj(ff,cond).data.nfr;
groupAv(cond).time = dataSbj(ff,cond).data.time;
groupAv(cond).pval = NaN(size(dataSbj(ff,cond).data.pval));
groupAv(cond).confradius = NaN(size(dataSbj(ff,cond).data.confradius));
groupAv(cond).i1f1 = dataSbj(ff,cond).data.i1f1;
groupAv(cond).elec = dataSbj(ff,cond).data.elec;
groupAv(cond).dtms = dataSbj(ff,cond).data.dtms;
groupAv(cond).nchan = dataSbj(ff,cond).data.nchan;
groupAv(cond).ndft = dataSbj(ff,cond).data.ndft;
groupAv(cond).freq = dataSbj(ff,cond).data.freq;
groupAv(cond).label = dataSbj(ff,cond).data.label;

% average
groupAv(cond).wave = mean(dataWave(:,:,cond,:),4);
groupAv(cond).sin = mean(dataSin(:,:,cond,:),4);
groupAv(cond).cos = mean(dataCos(:,:,cond,:),4);

groupAv(cond).amp = sqrt(groupAv(cond).sin.^2 + groupAv(cond).cos.^2);

groupAv(cond).filteredWave = mean(dataFil(:,:,cond,:),4);
groupAv(cond).noiseWave = mean(dataNoise(:,:,cond,:),4);

end


end
