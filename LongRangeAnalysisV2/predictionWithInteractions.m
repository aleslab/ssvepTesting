

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/AxxFilesAverageRef/';
% dataDir = '/Users/marleneponcet/Documents/dataLongRangeV2/AxxFiles/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    fprintf('Processing sbj %s \n',num2str(ff))
    
    % original LR motion
    sbj(ff,1).data = Axx.cond(1);
    
    % linear prediction
    sbj(ff,2).data = sumAxxWithShift(Axx.cond(3),Axx.cond(2)); %1st is the one to shift, 2nd does not 
    
    % linear prediction + spatial interaction
    % reconstructed LR simultaneous
    linearSpatial = sumAxx(Axx.cond(3),Axx.cond(2));   
    % original LR simultaneous
    actualSpatial = Axx.cond(4);  
    % difference
    nonLinearSpatial = computeDiff(actualSpatial,linearSpatial); % data1-data2
    % get half of the interaction term
    halfNonLinearSpatial = multiplyAxx(nonLinearSpatial,0.5);
    % spatial prediction
    sbj(ff,3).data = sumAxxWithShift(sumAxx(Axx.cond(3),halfNonLinearSpatial),sumAxx(Axx.cond(2),halfNonLinearSpatial));
    
    % linear prediction + temporal interaction
    % reconstructed LR flash left
    linearTemp1 = sumAxxWithShift(Axx.cond(2),Axx.cond(2));   
    % original LR flash left
    actualTemp1 = Axx.cond(5);  
    % difference
    nonLinearL = computeDiff(actualTemp1,linearTemp1);
     % get half of the interaction term
    halfNonLinearL = multiplyAxx(nonLinearL,0.5);
    % reconstructed LR flash right
    linearTemp2 = sumAxxWithShift(Axx.cond(3),Axx.cond(3));   
    % original LR flash right
    actualTemp2 = Axx.cond(6);  
    % difference and * 0.5
    halfNonLinearR = multiplyAxx(computeDiff(actualTemp2,linearTemp2),0.5);    
    % temporal prediction
    sbj(ff,4).data = sumAxxWithShift(sumAxx(Axx.cond(3),halfNonLinearR),sumAxx(Axx.cond(2),halfNonLinearL));
    
    % linear prediction + spatial interaction + temporal interaction 
    poolNonLinearL = multiplyAxx(sumAxx(halfNonLinearSpatial,nonLinearL),0.5);
    poolNonLinearR = multiplyAxx(sumAxx(halfNonLinearSpatial,halfNonLinearR),0.5);
    sbj(ff,5).data = sumAxxWithShift(sumAxx(Axx.cond(3),poolNonLinearR),sumAxx(Axx.cond(2),poolNonLinearL));
    
    % original SR motion
    sbj(ff,6).data = Axx.cond(7);
    
    % linear prediction
    sbj(ff,7).data = sumAxxWithShift(Axx.cond(9),Axx.cond(8)); %1st is the one to shift, 2nd does not 
    
    % linear prediction + spatial interaction
    % reconstructed SR simultaneous
    linearSpatial = sumAxx(Axx.cond(8),Axx.cond(9));   
    % original SR simultaneous
    actualSpatial = Axx.cond(10);  
    % difference
    nonLinearSpatial = computeDiff(actualSpatial,linearSpatial); % data1-data2
    % get half of the interaction term
    halfNonLinearSpatial = multiplyAxx(nonLinearSpatial,0.5);
    % spatial prediction
    sbj(ff,8).data = sumAxxWithShift(sumAxx(Axx.cond(9),halfNonLinearSpatial),sumAxx(Axx.cond(8),halfNonLinearSpatial));
    
    % linear prediction + temporal interaction
    % reconstructed SR flash one side (right?)
    linearTemp1 = sumAxxWithShift(Axx.cond(8),Axx.cond(8));   
    % original LR flash right
    actualTemp1 = Axx.cond(11);  
    % difference
    nonLinearL = computeDiff(actualTemp1,linearTemp1);
     % get half of the interaction term
    halfNonLinearL = multiplyAxx(nonLinearL,0.5);
    % reconstructed SR flash one side (left?)
    linearTemp2 = sumAxxWithShift(Axx.cond(9),Axx.cond(9));   
    % original LR flash left
    actualTemp2 = Axx.cond(12);  
    % difference and * 0.5
    halfNonLinearR = multiplyAxx(computeDiff(actualTemp2,linearTemp2),0.5);    
    % temporal prediciton
    sbj(ff,9).data = sumAxxWithShift(sumAxx(Axx.cond(9),halfNonLinearR),sumAxx(Axx.cond(8),halfNonLinearL));
    
    % linear prediction + spatial interaction + temporal interaction 
    poolNonLinearL = multiplyAxx(sumAxx(halfNonLinearSpatial,nonLinearL),0.5);
    poolNonLinearR = multiplyAxx(sumAxx(halfNonLinearSpatial,halfNonLinearR),0.5);
    sbj(ff,10).data = sumAxxWithShift(sumAxx(Axx.cond(9),poolNonLinearR),sumAxx(Axx.cond(8),poolNonLinearL));

    % on-linear spatio-temporal component
    % long-range
    sbj(ff,11).data = computeDiff(sbj(ff,1).data, sbj(ff,5).data);
    % short range
    sbj(ff,12).data = computeDiff(sbj(ff,6).data, sbj(ff,10).data);    
    
end

% do the average
gpCorrection = averageSbj(sbj);

%%%%
gpCorrection(1).condLabel = 'LR_originalmotion';
gpCorrection(2).condLabel = 'LR_linearPred';
gpCorrection(3).condLabel = 'LR_spatialPred';
gpCorrection(4).condLabel = 'LR_tempPred';
gpCorrection(5).condLabel = 'LR_spat&tempPred';
gpCorrection(6).condLabel = 'SR_originalmotion';
gpCorrection(7).condLabel = 'SR_linearPred';
gpCorrection(8).condLabel = 'SR_spatialPred';
gpCorrection(9).condLabel = 'SR_tempPred';
gpCorrection(10).condLabel = 'SR_spat&tempPred';
gpCorrection(11).condLabel = 'LR_SPnl';
gpCorrection(12).condLabel = 'SR_SPnl';

save('correctionAvRef','gpCorrection','cfg')





