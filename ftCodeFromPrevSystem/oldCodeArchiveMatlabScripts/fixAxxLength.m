function []=fixOddSpectrum(projectInfo,optns)
%function pdOdd2Eeglab(projectInfo)
%
outputDir = projectInfo.powerDivaExportDir;

allAxx = dir(fullfile(projectInfo.powerDivaExportDir,'Axx_*.mat'));

%Parse all raw files for number of conditions and trials
%Note tricky use of [] to turn struct into 1 long string;

if isempty([allAxx.name])
    error(['Cannot find RAW, single trial exports, for subject: ' projectInfo.subjId])
end

condAndTrialNum = sscanf([allAxx.name],'Axx_c%d.mat');

%Finds unique condition numbers in the export directory.
%condNum = unique(condAndTrialNum(1:2:end))';

condNum = unique(condAndTrialNum)'

%responseData = cell(length(condNum),1);
responseData = struct([]);

for iCond = condNum,
   
    iCond
    
    %[dat respTime condInfo] = sortOddStep(projectInfo.powerDivaExportDir,optns,iCond);

    
    filename = sprintf('Axx_c%.3i.mat',iCond);
    
    fileList = dir(fullfile(outputDir,filename));
    
    if length(fileList)==0,
        
        error(['No file found for condition: ' num2str(iCond) ] )
    end
    
    
    axxFileName = fullfile(outputDir,filename)
    axx = load(axxFileName);
    
    axx.Wave = axx.Wave(1:691,:);
	axx.nT = 691;

       

    % Create frequency domain data
	axx.dTms = 2.3125;
    
    axx.nFr = fix(size(axx.Wave,1)/2);
    dft = dftmtx(axx.nT);   

    dftDat = dft*axx.Wave;
    dftDat = dftDat(1:axx.nFr,:);
    
    
	%WARNING: freq domain display won't be right
    axx.Amp = abs(dftDat/axx.nFr);
    axx.Cos = real(dftDat)/axx.nFr;
    axx.Sin = -imag(dftDat)/axx.nFr;
    
	
	

    save(axxFileName,'-struct', 'axx');

end


