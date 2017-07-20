function pdCPT2Axx(projectInfo,optns)
%function pdCPT2Axx(projectInfo)
%
%Function to take powerdiva raw data and read it into EEGlab.
%
%


% Axx fields to fake:
% 
%          cndNmb: 2
%            nTrl: 25
%             nCh: 128
%              nT: 432
%             nFr: 108
%            dTms: 2.3125
%            dFHz: 0.5000
%            i1F1: 2
%            i1F2: 72
%     DataUnitStr: 'microVolts'
%            Wave: [432x128 double]
%             Cos: [108x128 double]
%             Sin: [108x128 double]
%             Amp: [108x128 double]
%             Cov: [128x128 double]

%outputDir = fullfile(projectInfo.projectDir,projectInfo.subjId,'_dev_');
outputDir = projectInfo.powerDivaExportDir;


%
allRaw = dir(fullfile(projectInfo.powerDivaExportDir,'Raw_*.mat'));

%Parse all raw files for number of conditions and trials
%Note tricky use of [] to turn struct into 1 long string;
condAndTrialNum = sscanf([allRaw.name],'Raw_c%d_t%d.mat');

%Finds unique condition numbers in the export directory.
condNum = unique(condAndTrialNum(1:2:end))';


% Axx fields to fake:
% 
%          cndNmb: 2
%            nTrl: 25
%             nCh: 128
%              nT: 432
%             nFr: 108
%            dTms: 2.3125
%            dFHz: 0.5000
%            i1F1: 2
%            i1F2: 72
%     DataUnitStr: 'microVolts'
%            Wave: [432x128 double]
%             Cos: [108x128 double]
%             Sin: [108x128 double]
%             Amp: [108x128 double]
%             Cov: [128x128 double]

  

for iCond = condNum,

%     condFilename = sprintf('Raw_c%0.3d_t*.mat',iCond);
%   
%     rawList = dir(fullfile(projectInfo.powerDivaExportDir,condFilename));
%     
%     for iRaw = 1:length(rawList),
%         
%         rawFiles{iRaw} = fullfile(projectInfo.powerDivaExportDir,rawList(iRaw).name);
%     end
% 
    
    iCond
%    [data info] = loadPDraw(rawFiles,0);
    optns.lowpass = false;
	optns.lock2reaction = true;
    clear dat paddedDat shiftedDat;

    [dat respTime condInfo] = eegResponseTime(projectInfo.powerDivaExportDir,optns,iCond);
	
	rtTarget = 600; %sample to allign all trials to
	respArray = respTime{1};
	respShift = respTime{1} - rtTarget;
	allTrials = length(respArray);
	for i=1:allTrials,
		if isnan(respShift(i))
			respShift(i) = 0;
		end
	end
	%shift Array can be zeros or NaNs
	shiftArray = zeros(size(dat{1}));
	%shiftArray(:,:,:) = NaN;
	
	
	
	
	
	
	
	
	
	
	
  
%    goodTrials = condInfo{1}.goodTrialLabel;

    % Axx fields that are correct:
    %
    Axx.cndNmb      = iCond;
    Axx.nTrl        = size(dat{1}, 2);    
    Axx.nCh         = size(dat{1},3);
    Axx.nT          = size(dat{1},1);
    Axx.dTms        = 10^3*(condInfo{1}.FreqHz)^-1;

%     for iTrial = 1:size(goodTrials,1),        
%         %if channel for this trial is ~good set to NaN
%         dat{1}(:,iTrial,~goodTrials(iTrial,:)) = NaN;
%     end
%     		

    
    % Creat Axx.Wave data
    if optns.lock2reaction == true
		
		for iT = 1:allTrials,
			iShift = respShift(iT);
			if iShift >= 0
				shiftArray(1:864-iShift,iT,:) = dat{1}(iShift+1:864,iT,:);
			else
				shiftArray(abs(iShift)+1:864,iT,:) = dat{1}(1:864+iShift,iT,:);
			end
			iT = iT +1;
		end
	          
		shiftDat = nanmean(shiftArray,2);
		shiftDat = squeeze(shiftDat);
    
    
        Axx.Wave        = shiftDat;
    else
        
         
        Axx.Wave        = squeeze(nanmean(dat{1},2));
   
    end
    
        
    %Hacky multiplier, fix this ---JMA
    Axx.Wave = 1e6*Axx.Wave;
    Axx.DataUnitStr = 'microVolts';

    % Create frequency domain data
    
    Axx.nFr = size(Axx.Wave,1);
    dft = dftmtx(Axx.nFr);   

    dftDat = dft*Axx.Wave;
    dftDat = dftDat(1:Axx.nFr,:);
    
    Axx.dFHz = (condInfo{1}.FreqHz)/size(Axx.Wave,1);

    % i1f1 = 3, means 1 Hz because of 3 oddsteps
    % Stuff tends to be multiples of 1 Hz. But we don't have that info
    % here so we are just going to set the index to 1 to make nF1 be all
    % freqs
    Axx.i1F1       = 1;
    Axx.i1F2       = 0;

    Axx.Amp = abs(dftDat);
    Axx.Cos = real(dftDat);
    Axx.Sin = imag(dftDat);
    
    %This is faked:
    Axx.Cov = eye(Axx.nCh);

    
      if optns.lock2reaction == true

          condNum2Write = 200+iCond;
      else
          condNum2Write = 100+iCond;
      end
         
      Axx.cndNmb = condNum2Write;
      
      filename=['Axx_c' num2str(condNum2Write,'%0.3d') '.mat'];
      file2write = fullfile(projectInfo.powerDivaExportDir,filename);
      
      
      save(file2write,'-struct','Axx')
    


end


