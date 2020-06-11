
function simulateActivityECC(projDir)

% This script sets up simulation paramaters and calls mrSimScript.

%The full path to the project. 
%Needs to be changed for unix/mac/windows
% projDir = '/Users/marleneponcet/Documents/data/topoMapForBiosemi';


%% Parameters using constraint model: left activation of dorsal+ventral (pooled)
%This gets the default set of simulation paramaters
params = skeriDefaultSimParameters;

% roiList = getRoisByType(roiDir,'func')
params.activeRoiList = {'V1Vecc-L','V1Decc-L','V2Vecc-L','V2Decc-L',...
    'V3Vecc-L','V3Decc-L','V3A-L','V4-L',...
    'MT-L','IPS-L','LOC-L'};  

%Left activation
%1f1 -> V1D+V1V
%2f1 -> V2D + V2V
%3f1 -> V3D+V3V
%4f1 -> V3A
%5f1 -> V4
%6f1 -> MT
%7f1 -> IPS
%8f1 -> LOC
params.roiHarm = {1e-4 1e-4 [0 1e-4] [0 1e-4]...
    [0 0 1e-4] [0 0 1e-4] [0 0 0 1e-4] [0 0 0 0 1e-4] ...
    [0 0 0 0 0 1e-4] [0 0 0 0 0 0 1e-4] [0 0 0 0 0 0 0 1e-4]};
params.stepTimeByRoi = true;
    
%Set this for each     
params.condNumber = 502;
mrSimScript(projDir,params);


%% Parameters using constraint model: left activation without pooling dorsal and ventral
%This gets the default set of simulation paramaters
params = skeriDefaultSimParameters;

% roiList = getRoisByType(roiDir,'func')
params.activeRoiList = {'V1Vecc-L','V1Decc-L','V2Vecc-L','V2Decc-L',...
    'V3Vecc-L','V3Decc-L','V3A-L','V4-L',...
    'MT-L','IPS-L','LOC-L'};  

%Left activation
%1f1 -> V1D
%2f1 -> V1V
%3f1 -> V2D 
%4f1 -> V2V
%5f1 -> V3D
%6f1 -> V3V
%7f1 -> V3A
%8f1 -> V4
%9f1 -> MT
%10f1 -> IPS
%11f1 -> LOC
params.roiHarm = {1e-4 [0 1e-4] [0 0 1e-4] [0 0 0 1e-4]...
       [0 0 0 0 1e-4] [0 0 0 0 0 1e-4] [0 0 0 0 0 0 1e-4] ...
       [0 0 0 0 0 0 0 1e-4] [0 0 0 0 0 0 0 0 1e-4] ...
       [0 0 0 0 0 0 0 0 0 1e-4] [0 0 0 0 0 0 0 0 0 0 1e-4]};
params.stepTimeByRoi = true;
    
%Set this for each     
params.condNumber = 503;
mrSimScript(projDir,params);


