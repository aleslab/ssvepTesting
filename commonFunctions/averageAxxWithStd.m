function [groupAv,sbjProj] = averageAxxWithStd(dataSbj)
% groupAv contains wave, sin, cos, amp average per condition
% sbjProj contains projected freq amplitudes for each participant per condtion (can be used to calculated std of amplitudes) 

for cond=1:size(dataSbj,1) 
% values that do not change
groupAv(cond).nt = dataSbj(cond,1).nt;
groupAv(cond).nfr = dataSbj(cond,1).nfr;
groupAv(cond).time = dataSbj(cond,1).time;
groupAv(cond).pval = NaN(size(dataSbj(cond,1).pval));
groupAv(cond).confradius = NaN(size(dataSbj(cond,1).confradius));
groupAv(cond).i1f1 = dataSbj(cond,1).i1f1;
groupAv(cond).elec = dataSbj(cond,1).elec;
groupAv(cond).dtms = dataSbj(cond,1).dtms;
groupAv(cond).nchan = dataSbj(cond,1).nchan;
groupAv(cond).ndft = dataSbj(cond,1).ndft;
groupAv(cond).freq = dataSbj(cond,1).freq;
groupAv(cond).label = dataSbj(cond,1).label;
groupAv(cond).condLabel = dataSbj(cond,1).condLabel;
end

% average
for cond=1:size(dataSbj,1) 
    clear dataTmpWave dataTmpSin dataTmpCos
    for ss=1:size(dataSbj,2)
        dataTmpWave(:,:,ss) = [dataSbj(cond,ss).wave];
        dataTmpSin(:,:,ss) = [dataSbj(cond,ss).sin];
        dataTmpCos(:,:,ss) = [dataSbj(cond,ss).cos];
        dataTmpAmp(:,:,ss) = [dataSbj(cond,ss).amp]; % sbj amplitude
        dataTmpPhi(:,:,ss) = atan2(dataSbj(cond,ss).cos,dataSbj(cond,ss).sin); % sbj phase
    end
    groupAv(cond).wave = mean(dataTmpWave,3);
    groupAv(cond).sin = mean(dataTmpSin,3);
    groupAv(cond).cos = mean(dataTmpCos,3);
    groupAv(cond).amp = sqrt(groupAv(cond).sin.^2 + groupAv(cond).cos.^2); % this is just to check that mean in group is the same as mean in the sbjProj
    groupTMPtheta = atan2(groupAv(cond).cos,groupAv(cond).sin); % phase of mean
    
    % projected amplitude for each sbj so that sdt can be computed
    % = individual subject amplitude projected onto mean vector
    % = individual participant amplitude scales by the cosine of the angle between the individual participant response and the group average
    for ss=1:size(dataSbj,2)
        sbjProj(:,:,cond,ss) = dataTmpAmp(:,:,ss) .* cos( dataTmpPhi(:,:,ss) - groupTMPtheta(:,:)); 
    end
    
end


end


