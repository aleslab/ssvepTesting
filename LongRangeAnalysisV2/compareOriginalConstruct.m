

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% gpData from 1 to 16:
% Long-range 
% reconstructed LR motion
% original LR motion
% reconstructed LR simultaneous
% original LR simultaneous
% reconstructed LR flash left
% original LR flash left
% reconstructed LR flash right
% original LR flash right

% short-range 
% reconstructed SR motion
% original SR motion
% reconstructed SR simultaneous
% original SR simultaneous
% reconstructed SR flash left
% original SR flash left
% reconstructed SR flash right
% original SR flash right

%%%%%%%%% gpDataDiff from 1 to 8:
% get differences between all the reconstructed and original
% LR motion
% LR simultaneous
% LR flash left
% LR flash right
% SR motion
% SR simultaneous
% SR flash left
% SR flash right


clear all;
addpath /Users/marleneponcet/Documents/Git/fieldtrip20170924
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
save('gpCompare','gpData','gpDataDiff','cfg')





