function [RT] = pdRawRT(projectInfo)
%function pdRawRT(projectInfo)
%spits out mean response time for subject
%Function to take powerdiva raw data and read it into EEGlab.
%
%
%for SKERI project
subjectList = dir(fullfile(projectInfo.projectDir, 'skeri*'));
%for Headless project
if isempty(subjectList)
subjectList = dir(fullfile(projectInfo.projectDir, 'headless*'));
end


for n = 1:length(subjectList)
	RT{n}.name = 'foo';
	RT{n}.responseTime = 0;
	RT{n}.falseAlarms = 0;
end

for iSubj = 1:length(subjectList),
    
	projectDir = projectInfo.projectDir;
    subjId = subjectList(iSubj).name;
    disp(['Processing subject: ' subjId ])
	PDname = dir(fullfile(projectDir,subjId,'Exp_MATL_*'));
    powerDivaExportDir = fullfile(projectDir,subjId,PDname(1).name);
    projectInfo.subjId = subjId;
    projectInfo.currentDir = fullfile(projectDir,subjId);
    projectInfo.powerDivaExportDir = powerDivaExportDir;
	
	allRaw = dir(fullfile(projectInfo.powerDivaExportDir,'Raw_*.mat'));

	%Parse all raw files for number of conditions and trials
	%Note tricky use of [] to turn struct into 1 long string;
	condAndTrialNum = sscanf([allRaw.name],'Raw_c%d_t%d.mat');

	%Finds unique condition numbers in the export directory.
% 	if isempty(condNum)
% 		condNum = unique(condAndTrialNum(1:2:end))';
% 	end
	
	condNum = unique(condAndTrialNum(1:2:end))';

	numCond = length(condNum);

	for iCond = 1:numCond,

		optns.lowpass = false;
		[respTime condInfo] = eegRespTime(projectInfo.powerDivaExportDir,optns,iCond);
		cycleLength = condInfo{1}.CycleLen;
		msLength = 1000*(condInfo{1}.CycleLen/condInfo{1}.FreqHz);
		if(sum(isnan(respTime{1}))>10)
			RT{iSubj}.falseAlarms = RT{iSubj}.falseAlarms + sum(~isnan(respTime{1}));
		else
		RT{iSubj}.responseTime(iCond, 1) = ((nanmean(respTime{1}))/cycleLength)* msLength;	
		end
		clear respTime condInfo;
		
	end
	clear respTime condInfo;
	
	RT{iSubj}.name = subjId;

end


