function [baselineDiff,baselineCI] = compute_baseline(dataWave,fixCond,compare)

clear baselineDiff baselineSTD baselineCI;
cc=0;
for cond = compare
    cc=cc+1; clear meanDiff;
    fprintf('compare with condition %d\n',cond)
    for tt=1:1000
        clear diffWave;
        for ss=1:size(dataWave,4)
            ord = Shuffle([fixCond cond]); 
            diffWave(:,:,ss) = dataWave(:,:,ord(1),ss) - dataWave(:,:,ord(2),ss);
        end
        meanDiff(tt,:,:) = mean(diffWave,3);
    end
    baselineDiff(:,:,cc) = squeeze(mean(meanDiff));
    baselineSTD(:,:,cc) = squeeze(std(meanDiff));
    baselineCI(:,:,cc) = 1.96 * (baselineSTD(:,:,cc) / sqrt(size(dataWave,4)));
end
end