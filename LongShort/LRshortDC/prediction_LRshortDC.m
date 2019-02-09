

clearvars;
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults

% dataDir = '/Users/marleneponcet/Documents/data/LRshortDC/V1/Axx/';
dataDir = '/Users/marleneponcet/Documents/data/LRshortDC/V2/Axx/';
outDir = '/Users/marleneponcet/Documents/data/LRshortDC/V2/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

for ff=1:length(listData)
    clear Axx;
    load([dataDir listData(ff).name]);
    fprintf('Processing sbj %s \n',num2str(ff))
    
    % original LR motion
    sbj(ff,1).data = Axx(1);
    
    % linear prediction
    sbj(ff,2).data = sumAxxWithShift(Axx(2),Axx(3)); %1st is the one to shift, 2nd does not 
    % difference recorded motion - linear prediction
    sbjDiff(ff,1).data = computeDiff(sbj(ff,1).data,sbj(ff,2).data); 
    
    % linear prediction + spatial interaction
    % reconstructed LR simultaneous
    linearSpatial = sumAxx(Axx(3),Axx(2));   
    % original LR simultaneous
    actualSpatial = Axx(4);  
    % difference
    nonLinearSpatial = computeDiff(actualSpatial,linearSpatial); % data1-data2
    interaction(ff,1).data = nonLinearSpatial;
    % get half of the interaction term
    halfNonLinearSpatial = multiplyAxx(nonLinearSpatial,0.5);
    % spatial prediction
    sbj(ff,3).data = sumAxxWithShift(sumAxx(Axx(2),halfNonLinearSpatial),sumAxx(Axx(3),halfNonLinearSpatial));
    % difference recorded motion - linear+spatial prediction
    sbjDiff(ff,2).data = computeDiff(sbj(ff,1).data,sbj(ff,3).data); 
    
    % linear prediction + temporal interaction
    % reconstructed LR flash left
    linearTemp1 = sumAxxWithShift(Axx(2),Axx(2));   
    % original LR flash left
    actualTemp1 = Axx(5);  
    % difference
    nonLinearL = computeDiff(actualTemp1,linearTemp1);
    interaction(ff,2).data = nonLinearL;
     % get half of the interaction term
    halfNonLinearL = multiplyAxx(nonLinearL,0.5);
    % reconstructed LR flash right
    linearTemp2 = sumAxxWithShift(Axx(3),Axx(3));   
    % original LR flash right
    actualTemp2 = Axx(6);  
    % difference and * 0.5
    interaction(ff,3).data = computeDiff(actualTemp2,linearTemp2);
    halfNonLinearR = multiplyAxx(computeDiff(actualTemp2,linearTemp2),0.5);    
    % temporal prediction
    sbj(ff,4).data = sumAxxWithShift(sumAxx(Axx(2),halfNonLinearL),sumAxx(Axx(3),halfNonLinearR));
    % difference recorded motion - linear+temporal prediction
    sbjDiff(ff,3).data = computeDiff(sbj(ff,1).data,sbj(ff,4).data); 

    % linear prediction + spatial interaction + temporal interaction 
    poolNonLinearL_LR = sumAxx(halfNonLinearSpatial,halfNonLinearL);
    poolNonLinearR_LR = sumAxx(halfNonLinearSpatial,halfNonLinearR);
    sbj(ff,5).data = sumAxxWithShift(sumAxx(Axx(2),poolNonLinearL_LR),sumAxx(Axx(3),poolNonLinearR_LR));
     % difference recorded motion - full prediction
    sbjDiff(ff,4).data = computeDiff(sbj(ff,1).data,sbj(ff,5).data); 
    
    % original SR motion
    sbj(ff,6).data = Axx(7);
    
    % linear prediction
    sbj(ff,7).data = sumAxxWithShift(Axx(8),Axx(9)); %1st is the one to shift, 2nd does not 
    % difference recorded motion - prediction
    sbjDiff(ff,5).data = computeDiff(sbj(ff,6).data,sbj(ff,7).data); 
    
    % linear prediction + spatial interaction
    % reconstructed SR simultaneous
    linearSpatial = sumAxx(Axx(8),Axx(9));   
    % original SR simultaneous
    actualSpatial = Axx(10);  
    % difference
    nonLinearSpatial = computeDiff(actualSpatial,linearSpatial); % data1-data2
    interaction(ff,4).data = nonLinearSpatial;
    % get half of the interaction term
    halfNonLinearSpatial = multiplyAxx(nonLinearSpatial,0.5);
    % spatial prediction
    sbj(ff,8).data = sumAxxWithShift(sumAxx(Axx(8),halfNonLinearSpatial),sumAxx(Axx(9),halfNonLinearSpatial));
    % difference recorded motion - prediction
    sbjDiff(ff,6).data = computeDiff(sbj(ff,6).data,sbj(ff,8).data); 
    
    % linear prediction + temporal interaction
    % reconstructed SR flash one side (right?)
    linearTemp1 = sumAxxWithShift(Axx(8),Axx(8));   
    % original LR flash right
    actualTemp1 = Axx(11);  
    % difference
    nonLinearL = computeDiff(actualTemp1,linearTemp1);
    interaction(ff,5).data = nonLinearL;
     % get half of the interaction term
    halfNonLinearL = multiplyAxx(nonLinearL,0.5);
    % reconstructed SR flash one side (left?)
    linearTemp2 = sumAxxWithShift(Axx(9),Axx(9));   
    % original LR flash left
    actualTemp2 = Axx(12);  
    % difference and * 0.5
    interaction(ff,6).data = computeDiff(actualTemp2,linearTemp2);
    halfNonLinearR = multiplyAxx(computeDiff(actualTemp2,linearTemp2),0.5);    
    % temporal prediction
    sbj(ff,9).data = sumAxxWithShift(sumAxx(Axx(8),halfNonLinearL),sumAxx(Axx(9),halfNonLinearR));
    % difference recorded motion - prediction
    sbjDiff(ff,7).data = computeDiff(sbj(ff,6).data,sbj(ff,9).data); 
    
    % linear prediction + spatial interaction + temporal interaction 
    poolNonLinearL_SR = sumAxx(halfNonLinearSpatial,halfNonLinearL);
    poolNonLinearR_SR = sumAxx(halfNonLinearSpatial,halfNonLinearR);
    sbj(ff,10).data = sumAxxWithShift(sumAxx(Axx(8),poolNonLinearL_SR),sumAxx(Axx(9),poolNonLinearR_SR));
    % difference recorded motion - prediction
    sbjDiff(ff,8).data = computeDiff(sbj(ff,6).data,sbj(ff,10).data); 
    
    % Not "real" non-linear component but useful for checks 
    interaction(ff,7).data = poolNonLinearL_LR;
    interaction(ff,8).data = poolNonLinearR_LR;
    interaction(ff,9).data = poolNonLinearL_SR;
    interaction(ff,10).data = poolNonLinearR_SR;    

    % non-linear spatio-temporal component
    % long-range
    sbj(ff,11).data = computeDiff(sbj(ff,1).data, sbj(ff,5).data);
    % short range
    sbj(ff,12).data = computeDiff(sbj(ff,6).data, sbj(ff,10).data);    
end


save([outDir 'sbjprediction'],'sbj','cfg')
save([outDir 'sbjDiff'],'sbjDiff','cfg')
save([outDir 'NLinteraction'],'interaction','cfg')



figure;hold on;
plot(halfNonLinearL.wave(23,:),'LineWidth',2);
plot(halfNonLinearR.wave(23,:),'LineWidth',2);
plot(halfNonLinearSpatial.wave(23,:),'LineWidth',2);
plot(poolNonLinearL_SR.wave(23,:),'LineWidth',2);
plot(poolNonLinearR_SR.wave(23,:),'LineWidth',2);
ylim([-6 6])
legend('0.5Tl','0.5Tr','0.5S','suml','sumr')

figure;hold on;
plot(nonLinearL.wave(23,:),'LineWidth',2);
plot(interaction(13,6).data.wave(23,:),'LineWidth',2);
plot(nonLinearSpatial.wave(23,:),'LineWidth',2);
plot(poolNonLinearL_SR.wave(23,:),'LineWidth',2);
plot(poolNonLinearR_SR.wave(23,:),'LineWidth',2);
ylim([-6 6])
legend('Tl','Tr','S','suml','sumr')



% 
% % do the average
% gpPred = averageSbj(sbj);
% gpInteraction = averageAxx(permute(interaction,[2 1]));
% 
% %%%%
% gpPred(1).condLabel = 'LR_originalmotion';
% gpPred(2).condLabel = 'LR_linearPred';
% gpPred(3).condLabel = 'LR_spatialPred';
% gpPred(4).condLabel = 'LR_tempPred';
% gpPred(5).condLabel = 'LR_spat&tempPred';
% gpPred(6).condLabel = 'SR_originalmotion';
% gpPred(7).condLabel = 'SR_linearPred';
% gpPred(8).condLabel = 'SR_spatialPred';
% gpPred(9).condLabel = 'SR_tempPred';
% gpPred(10).condLabel = 'SR_spat&tempPred';
% % gpCorrection(11).condLabel = 'LR_SPnl';
% % gpCorrection(12).condLabel = 'SR_SPnl';
% 
% gpInteraction(1).condLabel = 'LRspatial';
% gpInteraction(2).condLabel = 'LRl';
% gpInteraction(3).condLabel = 'LRr';
% gpInteraction(7).condLabel = 'LRst';
% gpInteraction(4).condLabel = 'SRspatial';
% gpInteraction(5).condLabel = 'SRl';
% gpInteraction(6).condLabel = 'SRr';
% gpInteraction(8).condLabel = 'SRst';
% 
% 
% 
% 
% interactiveSteadyStatePlot2(cfg,gpPred)

