

load('correctionAvRef.mat')
load('sbjCorrectionAvRef.mat')

for ff=1:size(sbj,1)
    for cond=1:size(sbj,2)
        dataWave(:,:,cond,ff) = sbj(ff,cond).data.wave;
    end
end
% dataWave(elec,time,condition,subj)


for cond = 2:5
    fprintf('compare with condition %d\n',cond)
    for tt=1:1000
        clear diffWave;
        for sbjNb=1:16
            ord = Shuffle([1 cond]); ordSR = Shuffle([6 cond+5]);
            diffWave(:,:,sbjNb) = dataWave(:,:,ord(1),sbjNb) - dataWave(:,:,ord(2),sbjNb);
            diffWaveSR(:,:,sbjNb) = dataWave(:,:,ordSR(1),sbjNb) - dataWave(:,:,ordSR(2),sbjNb);
        end
        meanDiff(tt,:,:) = mean(diffWave,3);
        meanDiffSR(tt,:,:) = mean(diffWaveSR,3);
    end
    baselineDiff(:,:,cond-1) = squeeze(mean(meanDiff));
    baselineDiffSR(:,:,cond-1) = squeeze(mean(meanDiffSR));
    baselineSTD(:,:,cond-1) = squeeze(std(meanDiff));
    baselineSTDSR(:,:,cond-1) = squeeze(std(meanDiffSR));
end


condName = {'L','LS','LT','LST'};

elec = 19;elec=23;
% long range
figure; hold on;
for compCond = 1:4
    subplot(2,2,compCond); hold on;
    plot(gpCorrection(1).time,baselineDiff(elec,:,compCond),'k')
%         errorbar(gpCorrection(1).time,baselineDiff(elec,:,compCond),baselineSTD(elec,:,compCond),'Color','k','linestyle','none');
    fill([gpCorrection(1).time gpCorrection(1).time(end:-1:1) gpCorrection(1).time(1)],...
        [baselineDiff(elec,:,compCond)-baselineSTD(elec,:,compCond) ...
        baselineDiff(elec,end:-1:1,compCond)+baselineSTD(elec,end:-1:1,compCond) ...
        baselineDiff(elec,1,compCond)-baselineSTD(elec,1,compCond)], 'k','EdgeColor','none','facealpha',0.1);
    plot(gpCorrection(1).time,gpCorrection(1).wave(elec,:) - gpCorrection(compCond+1).wave(elec,:),'r')
    ylim([-1 1])
    title(condName{compCond})
end
% short range
figure; hold on;
for compCond = 1:4
    subplot(2,2,compCond); hold on;
    plot(gpCorrection(1).time,baselineDiffSR(elec,:,compCond),'k')
%     errorbar(gpCorrection(1).time,baselineDiffSR(elec,:,compCond+4),baselineSTDSR(elec,:,compCond+4),'Color','k','linestyle','none');
    fill([gpCorrection(1).time gpCorrection(1).time(end:-1:1) gpCorrection(1).time(1)],...
        [baselineDiffSR(elec,:,compCond)-baselineSTDSR(elec,:,compCond) ...
        baselineDiffSR(elec,end:-1:1,compCond)+baselineSTDSR(elec,end:-1:1,compCond) ...
        baselineDiffSR(elec,1,compCond)-baselineSTDSR(elec,1,compCond)], 'k','EdgeColor','none','facealpha',0.1);
    plot(gpCorrection(1).time,gpCorrection(6).wave(elec,:) - gpCorrection(compCond+6).wave(elec,:),'r')
    ylim([-1 1])
    title(condName{compCond})
end


