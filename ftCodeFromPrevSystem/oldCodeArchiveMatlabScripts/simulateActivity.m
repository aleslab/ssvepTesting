% This script sets up simulation paramaters and calls mrSimScript.

%The full path to the project. 
%Needs to be changed for unix/mac/windows
projDir = '/Volumes/MRI-1/4D2/JMA_PROJECTS/c1v1flip/mrcProj'
%projDir = '/raid/MRI/data/4D2/JMA_PROJECTS/c1v1flip/mrcProj'


%% Parameters for using V1/3/4 D/V
%This gets the default set of simulation paramaters
params = skeriDefaultSimParameters

%Use V1/2/3 D/V subdivisions
params.activeRoiList = { 'V1D-L' 'V1V-L' 'V1D-R' 'V1V-R' ...
    'V2D-L' 'V2V-L' 'V2D-R' 'V2V-R' ...
    'V3D-L' 'V3V-L' 'V3D-R' 'V3V-R'};  


%Bilat activation
%1f1 -> V1D
%2f1 -> V1V
%3f1 -> V2D
%4f1 -> V2V
%5f1 -> V3D
%6f1 -> V3V
params.roiHarm = { 1e-4 [0 1e-4] 1e-4 [0 1e-4] ...
    [0 0 1e-4] [0 0 0 1e-4] [0 0 1e-4] [0 0 0 1e-4] ...
    [0 0 0 0 1e-4] [0 0 0 0 0 1e-4] [0 0 0 0 1e-4] [0 0 0 0 0 1e-4]};
    
%Set this for each     
params.condNumber = 901;
mrSimScript(projDir,params);

%% Parameters for simulating both V1 and V2 waveform with 90 degree phase
%This gets the default set of simulation paramaters
%params = skeriDefaultSimParameters
params = jakeDefaultSimParameters

params.stepTimeByRoi = false;

% Bilat V1+V2 on 1f1
params.activeRoiList = {'V1-L'  'V1-R'  'V2D-L'  'V2D-R' 'V2V-L' 'V2V-R'};
params.roiHarm = { [1e-4] [1e-4] [complex(0,1e-4)] [complex(0,1e-4)] [complex(0,1e-4)] [complex(0,1e-4)] } 
params.condNumber = 902;

mrSimScript(projDir,params);

%% Parameters for simulating arbitrary time functions.
%This one is tricky. No error checking is done for the correct number of
%time points. So you must choose those carefully, and they
%MUST be consistent with the other Axx_ datasets!
%params = skeriDefaultSimParameters
params = jakeDefaultSimParameters

params.stepTimeByRoi = false;

%params.activeRoiList = {'V1-L'  'V1-R'  'V2D-L'  'V2D-R' 'V2V-L' 'V2V-R'};

params.activeRoiList = {'V1D-L'  'V1D-R'  'V1V-L'  'V1V-R' 'V2D-L'  'V2D-R' 'V2V-L' 'V2V-R'};

%roiHarm must be set to something. so I used the values from above
params.roiHarm = { [1e-4] [1e-4] [1e-4] [1e-4] [complex(0,1e-4)] [complex(0,1e-4)] [complex(0,1e-4)] [complex(0,1e-4)] } 


% %This initializes to wave forms for each source to 0;
for i=1:length(params.activeRoiList),
    params.roiTime{i} = zeros(100,1);
end

%Turn each source on for 10 time samples.
%1e-4 sources give voltages around ~microvolt scalp potentials
startIdx = 15;
timeOn = 10;


for iArea = 1:length(params.activeRoiList),
    
    
    thisStartIdx = startIdx + (iArea-1)*timeOn;
    theseTimes = [1:timeOn] + thisStartIdx;
    
    params.roiTime{iArea}(theseTimes) = 1e-4;

end


params.condNumber = 903;

params

%mrSimScript(projDir,params);

%% Simulation of just V1D
%params = skeriDefaultSimParameters
params = jakeDefaultSimParameters

params.stepTimeByRoi = true;

params.activeRoiList = {'V1D-L' 'V1D-R'};

params.roiHarm = {1e-4 1e-4} 

% for i=1:length(params.activeRoiList),
%     params.roiTime{i} = zeros(100,1);
% end
% 
% params.roiTime{1}(31:40) = 1e-4;
% params.roiTime{2}(41:50) = 1e-4;

params.condNumber = 904;

params

mrSimScript(projDir, params);

%% Simulation of V1V

%params = skeriDefaultSimParameters
params = jakeDefaultSimParameters

params.stepTimeByRoi = true;

params.activeRoiList = {'V1V-L' 'V1V-R'};

params.roiHarm = {1e-4 1e-4} 

params.condNumber = 905;

params

mrSimScript(projDir, params);


%% Simulation of V1D and V1V 

%params = skeriDefaultSimParameters
params = jakeDefaultSimParameters

params.stepTimeByRoi = true;

params.activeRoiList = {'V1D-L' 'V1D-R' 'V1V-L' 'V1V-R'};

params.roiHarm = {1e-4 1e-4 1e-4 1e-4} 

params.condNumber = 906;

params

mrSimScript(projDir, params);

%% Parameters for looping over v1/2/3
%This gets the default set of simulation paramaters
params = skeriDefaultSimParameters

%Use V1/2/3 D/V subdivisions
rois2Use = { ...
    'V1D-L' 'V1V-L' 'V2D-L' 'V2V-L' 'V3D-L' 'V3V-L' ...    
    'V1D-R' 'V1V-R' 'V2D-R' 'V2V-R' 'V3D-R' 'V3V-R'};

baseCndNumber = [951:956 961:966];




params.roiHarm = { 1e-4 };
    
for iRoi = 1:length(rois2Use),

    params.activeRoiList = rois2Use(iRoi);
    params.condNumber = baseCndNumber(iRoi);

    mrSimScript(projDir,params);
end


%% Comments
%The following are the default paramaters along with a description
%
% params.activeRoiList = {'V1-L'  'V1-R'  'LOC-L'  'LOC-R'};
% activeRoiList is a cell array of the names of the ROIs to simulate
% the ROI must exist in the subjects directory.
% 
% params.roiHarm = {[1]  [0 0 + 1.0000i]  [0 0 1]  [0 0 0 1]};
% roiHarm is a slightly tricky one. It is a cell array containing vectors
% of complex values that index into harmonics. The value is the current
% source density put in the roi
% e.g. [1]      -> 1 at 1f1 @ 0 phase angle
%      [1i]     -> 1 at 1f1 @ 90 phase angle
%      [ 0 1]   -> 1 at 2f1 @ 0 phase angle
%      [ 0 0 1i] -> 1 at 3f1 @ 90 phase angle
%
%
% params.condNumber = 901;
% condNumber is the name given to the axx file: Axx_c901.mat
%
%
% params.stepTimeByRoi = true;
% stepTimeByRoi makes a synthetic time function.
% If it is false the time domain data is constructed from the chosen
% harmonics
% if it is true the time domain is just uniform activation of each roi in
% order. 
%
% params.roiTime is not set by default
% Please read the following and be carefull with this parameter.
% if params.stepTimeByRoi is false, and this is non empty the simulator
% will use the waveforms in this cell array to simulate data
% This one is tricky. No error checking is done for the correct number of
% time points. So you must choose those carefully, and they
% MUST be consistent with the other Axx_ datasets!
% 
%
% params.sphereModel = false;
% sphereModel chooses which forward model to use.
% sphereModel = false means use a BEM model
% sphereModel = true means use a spherical head model
% caution: check for the file skeri????/_MNE_/skeri????-sph-fwd.fif
% to ensure spherical model calculations have been done.
% if not run: prepareInversesForMrc again
% 
%
%
% params.noise.type = 'white' 
% noise.type can be either 'white' or 'colored';
% white uses an identity covariance matrix
% colored uses the covariance estimate from powerDiva
%
% params.noise.level = 0;
% noise.level is calculated as a percentage the mean activation in time
% e.g. noise = noise.level*mean(abs(wave(:)))*randn(nT,nElec)




