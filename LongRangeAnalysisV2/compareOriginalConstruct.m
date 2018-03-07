
clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/AxxFiles/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};



for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    fprintf('Processing sbj %s \n',num2str(ff))
    
    %%%% LONG RANGE
    % reconstructed LR motion
    sbj(ff,1).data = reconstructTiming(Axx.cond(3),Axx.cond(2)); %1st is the one to shift, 2nd does not change
    % original LR motion
    sbj(ff,2).data  = Axx.cond(1);
    % difference
    sbjDiff(ff,1).data = computeDiff(sbj(ff,1).data,sbj(ff,2).data); % compute the difference (reconstructed always first)
    
    % reconstructed LR simultaneous
    sbj(ff,3).data = reconstructSimult(Axx.cond(3),Axx.cond(2));   
    % original LR simultaneous
    sbj(ff,4).data = Axx.cond(4);  
    % difference
    sbjDiff(ff,2).data = computeDiff(sbj(ff,3).data,sbj(ff,4).data); % compute the difference (reconstructed always first)
    
    % reconstructed LR flash left
    sbj(ff,5).data = reconstructTiming(Axx.cond(2),Axx.cond(2)); %1st is the one to shift, 2nd does not change    
    % original LR flash left
    sbj(ff,6).data = Axx.cond(5);  
    % difference
    sbjDiff(ff,3).data = computeDiff(sbj(ff,5).data,sbj(ff,6).data); % compute the difference (reconstructed always first)
    
    % reconstructed LR flash right
    sbj(ff,7).data = reconstructTiming(Axx.cond(3),Axx.cond(3)); %1st is the one to shift, 2nd does not change    
    % original LR flash right
    sbj(ff,8).data = Axx.cond(6);  
    % difference 
    sbjDiff(ff,4).data = computeDiff(sbj(ff,7).data,sbj(ff,8).data); % compute the difference (reconstructed always first)

    %%%% SHORT RANGE
    % reconstructed SR motion
    sbj(ff,9).data = reconstructTiming(Axx.cond(9),Axx.cond(8)); %1st is the one to shift, 2nd does not change
    % original SR motion
    sbj(ff,10).data  = Axx.cond(7);
    % difference
    sbjDiff(ff,5).data = computeDiff(sbj(ff,9).data,sbj(ff,10).data); % compute the difference (reconstructed always first)
    
    % reconstructed SR simultaneous
    sbj(ff,11).data = reconstructSimult(Axx.cond(9),Axx.cond(8));   
    % original SR simultaneous
    sbj(ff,12).data = Axx.cond(10);  
    % difference
    sbjDiff(ff,6).data = computeDiff(sbj(ff,11).data,sbj(ff,12).data); % compute the difference (reconstructed always first)

    % reconstructed LR flash left
    sbj(ff,13).data = reconstructTiming(Axx.cond(8),Axx.cond(8)); %1st is the one to shift, 2nd does not change    
    % original LR flash left
    sbj(ff,14).data = Axx.cond(11);  
    % difference
    sbjDiff(ff,7).data = computeDiff(sbj(ff,13).data,sbj(ff,14).data); % compute the difference (reconstructed always first)
    
    % reconstructed LR flash right
    sbj(ff,15).data = reconstructTiming(Axx.cond(9),Axx.cond(9)); %1st is the one to shift, 2nd does not change    
    % original LR flash right
    sbj(ff,16).data = Axx.cond(12);  
    % difference 
    sbjDiff(ff,8).data = computeDiff(sbj(ff,15).data,sbj(ff,16).data); % compute the difference (reconstructed always first)
end

% do the average
gpData = averageSbj(sbj);
gpDataDiff = averageSbj(sbjDiff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% gpData from 1 to 16:
% Long-range 
gpData(1).condLabel = 'reconstruct_LRmotion';
gpData(2).condLabel = 'original_LRmotion';
gpData(3).condLabel = 'reconstruct_LRsimult';
gpData(4).condLabel = 'original_LRsimult';
gpData(5).condLabel = 'reconstruct_LRflashLeft'; % probably right side
gpData(6).condLabel = 'original_LRflashLeft';
gpData(7).condLabel = 'reconstruct_LRflashRight';
gpData(8).condLabel = 'original_LRflashRight';

% short-range 
gpData(9).condLabel = 'reconstruct_SRmotion';
gpData(10).condLabel = 'original_SRmotion';
gpData(11).condLabel = 'reconstruct_SRsimult';
gpData(12).condLabel = 'original_SRsimult';
gpData(13).condLabel = 'reconstruct_SRflashLeft'; % probably right side
gpData(14).condLabel = 'original_SRflashLeft';
gpData(15).condLabel = 'reconstruct_SRflashRight';
gpData(16).condLabel = 'original_SRflashRight';

%%%%%%%%% gpDataDiff from 1 to 8:
% get differences between all the reconstructed and original
gpDataDiff(1).condLabel = 'LRmotion';
gpDataDiff(2).condLabel = 'LRsimult';
gpDataDiff(3).condLabel = 'LRflashLeft';
gpDataDiff(4).condLabel = 'LRflashRigth';
gpDataDiff(5).condLabel = 'SRmotion'; % probably right side
gpDataDiff(6).condLabel = 'SRsimult';
gpDataDiff(7).condLabel = 'SRflashLeft';
gpDataDiff(8).condLabel = 'SRflashRight';% LR motion

save('gpCompare','gpData','gpDataDiff','cfg')
save('sbjResults','sbjDiff','cfg')





