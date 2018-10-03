

singleResp1 = [0 1 1 0 1 0];
singleResp2 = [0 0 2 0 1 0];

temporalInteraction = [0 -.2 -.5 -.5 0 0];
spatialInteraction  = [0 -.3 0 -.1 -.1 0];

% Actual motion signal with only linear/spat/temp/spatTemp interactions
linear = singleResp1 + circshift(singleResp2,[0 3]);
spatAct = singleResp1 + spatialInteraction + circshift(singleResp2+spatialInteraction,[0 3]);
tempAct = (singleResp1 + temporalInteraction) + circshift(singleResp2+temporalInteraction,[0 3]);
spatTempAct = (singleResp1 + temporalInteraction + spatialInteraction) + circshift(singleResp2+temporalInteraction+ spatialInteraction,[0 3]);

% Predictions:

% Linear prediction
linearPred = singleResp1 + circshift(singleResp2,[0 3]);

% only spatial interaction
simultCondition = (singleResp1 +spatialInteraction) + (singleResp2 + spatialInteraction);
spatialLinearPred  = singleResp1 + singleResp2;
estimateSpatNonlinear = simultCondition-spatialLinearPred;
spatPred = singleResp1 + estimateSpatNonlinear/2 + circshift(singleResp2+estimateSpatNonlinear/2,[0 3]);

% only temporal interactions
temporalLinearPred1 = singleResp1 + circshift(singleResp1,[0 3]);
temporal1 = (singleResp1 + temporalInteraction)+circshift(singleResp1+temporalInteraction,[0 3]);
estimateTemporalNonlinear1 = temporal1-temporalLinearPred1;
 
temporalLinearPred2 = singleResp2 + circshift(singleResp2,[0 3]);
temporal2     = (singleResp2 + temporalInteraction)+circshift(singleResp2+temporalInteraction,[0 3]);
estimateTemporalNonlinear2 = temporal2-temporalLinearPred2;

tempPred = singleResp1 + estimateTemporalNonlinear1/2 + circshift(singleResp2+estimateTemporalNonlinear2/2,[0 3]);

% both spatial and temporal
spatTempPred = singleResp1 + estimateTemporalNonlinear1/2  + estimateSpatNonlinear/2 + circshift(singleResp2+estimateTemporalNonlinear2/2+ + estimateSpatNonlinear/2,[0 3]);




