function [groupAv,sbjProj] = averageAxxWithStd(dataSbj)
% groupAv contains wave, sin, cos, amp average per condition
% sbjProj contains projected freq amplitudes for each participant per condition (can be used to calculate std of amplitudes) 

% copy fields 
for cond=1:size(dataSbj,1)
    for fn = fieldnames(dataSbj)'
        groupAv(cond).(fn{1}) = dataSbj(cond,1).(fn{1});
    end
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


