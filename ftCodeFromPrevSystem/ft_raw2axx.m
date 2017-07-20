function [cfg] = ft_raw2axx(projectInfo,optns)
%function ft_raw2Axx(projectInfo)
%
%Function to take powerdiva raw data, use field trip to reject artifacts and export an Axx.
%

% Axx fields to fake:
% 
%          cndNmb: 2
%            nTrl: 70
%             nCh: 128
%              nT: 864
%             nFr: 108
%            dTms: 2.3125
%            dFHz: 0.5000
%            i1F1: 2
%            i1F2: 72
%     DataUnitStr: 'microVolts'
%            Wave: [864x128 double]
%             Cos: [108x128 double]
%             Sin: [108x128 double]
%             Amp: [108x128 double]
%             Cov: [128x128 double]


%projectInfo has the following fields
%projectInfo.projectDir -> the project directory currently being analyzed
%projectInfo.subjId -> the current subject, eg skeri0001
%projectInfo.currentDir -> full path to the current subjects directory
%projectInfo.powerDivaExportDir -> full path to directory containing PD
%exports;


inputDir = fullfile(projectInfo.projectDir,projectInfo.subjId,'Raw');
outputDir = projectInfo.powerDivaExportDir;


allRaw = dir(fullfile(inputDir,'*.raw'));
Raw = allRaw(length(allRaw)).name;

%currently needs elp file to run... it's not actually necessary for any of
%the artifact rejection
dataset = fullfile(inputDir, Raw);
elps = dir(fullfile(projectInfo.projectDir, projectInfo.subjId, 'Polhemus', '*.elp'));
polhemus = fullfile(projectInfo.projectDir, projectInfo.subjId, 'Polhemus', elps.name);


%read in a config and data structure for the EGI file
[cfg, data, condList] = ft_readEGI2(dataset, polhemus);

%particular to the CPT data 
%timing correction is the real sampling rate/ netstation's default sampling
%rate (1kHz)
warning('check sampling rate from powerDiva video')
timingCorrection = 432/1000;

%sort stucture by condition
for i=1:length(condList)
	tmpArray = strcmp(condList{i}, cfg.event);
	cfgConds{i} = cfg;
	cfgConds{i}.event = cfg.event(tmpArray);
	cfgConds{i}.trl = cfg.trl(tmpArray);
	cfgConds{i}.trialdef.eventvalue = condList{i};
	dataConds{i} = data;
	dataConds{i}.trial = data.trial(tmpArray);
	dataConds{i}.time = data.time(tmpArray);
	dataConds{i}.cfg = cfgConds{i};
end


for a = 1:length(condList)
	tmpWave = ft_timelockanalysis(cfgConds{a}, dataConds{a});
	%Axx{a}.wave = tmpWave;
	Axx.cndNmb = a;
    % Axx fields that are correct:
    %
    Axx.nTrl        = size(tmpWave.trial, 1);    
    Axx.nCh         = size(tmpWave.trial,2);
    Axx.nT          = size(tmpWave.trial,3);
    Axx.dTms        = 10^3*(timingCorrection*1000)^-1;		
    Axx.Wave = tmpWave.avg';
	
        
    %Hacky multiplier, fix this ---JMA
    %Axx.Wave = 1e6*Axx.Wave;
    Axx.DataUnitStr = 'microVolts';

    % Create frequency domain data
    
    Axx.nFr = fix(size(Axx.Wave,1)/2);
    dft = dftmtx(Axx.nT);   

    dftDat = dft*Axx.Wave;
    dftDat = dftDat(1:Axx.nFr,:);
    
    Axx.dFHz = (timingCorrection*1000)/size(Axx.Wave,1);

    % i1f1 = 3, means 1 Hz because of 3 oddsteps
    % Stuff tends to be multiples of 1 Hz. But we don't have that info
    % here so we are just going to set the index to 1 to make nF1 be all
    % freqs
    Axx.i1F1       = 1;
    Axx.i1F2       = 0;
	
	
	%WARNING: freq domain display won't be right
    Axx.Amp = abs(dftDat/Axx.nFr);
    Axx.Cos = real(dftDat)/Axx.nFr;
    Axx.Sin = -imag(dftDat)/Axx.nFr;
    
    %This is faked:
    Axx.Cov = eye(Axx.nCh);

    
    condNum2Write = a;
         
    Axx.cndNmb = condNum2Write;
      
      filename=['Axx_c' num2str(condNum2Write,'%0.3d') '.mat'];
      file2write = fullfile(projectInfo.powerDivaExportDir,filename);
      
      
      save(file2write,'-struct','Axx')
    


end



