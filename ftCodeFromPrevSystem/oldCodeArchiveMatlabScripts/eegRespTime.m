function [respTime condInfo] = eegRespTime(dataDir,optns,condition)
%function [dat respTime condInfo] = eegResponseTime(dataDir,optns)
%returns response times for button presses stored a channel in the EEG
%array

if nargin<=1 || isempty(optns)
    optns.lowpass = false;
end

if nargin<=2 || isempty(condition)
    cList = dir(fullfile(dataDir,'Raw_c*_t001.mat'));
    
    condAndTrialNum = sscanf([cList.name],'Raw_c%d_t%d.mat');

    %Finds unique condition numbers in the export directory.
    condList =  unique(condAndTrialNum(1:2:end))';
else
    condList = condition;

end


condIdx = 1;
for iCnd=condList,


    [condData] = loadPowerDivaRaw(dataDir,iCnd);
	respTimes = zeros(length(condData),1);
    idx = 1;
    clear fullData;
    idx =1;
	respTimes = zeros(length(condData),1);
	allTrials = length(condData);

    for iTrial = 1:length(condData),
		
		nEpochs = condData{iTrial}.NmbEpochs; % # epochs to extract
        nEeg    = condData{iTrial}.NmbChanEEG;% # of eeg channels to read
        nPrelude = condData{iTrial}.NmbPreludeEpochs;
		rTime = find(condData{iTrial}.RawTrial(:,131));
		if isempty(rTime)
			respTimes(iTrial) = NaN;
			%respTimes(iTrial) = 0;
		else
			respTimes(iTrial) = rTime(1,1);
		end
		
              

        if iTrial == 1
            nTotalEpochs = (nEpochs-nPrelude)*length(condData);
          
            CycleLen = size(condData{iTrial}.RawTrial,1)./condData{iTrial}.NmbEpochs;

            condInfo{condIdx} = rmfield(condData{iTrial},'RawTrial');
            condInfo{condIdx}.CycleLen = CycleLen;
		end

        validEpochs = 1:nEpochs;


        idx = idx+length(validEpochs);
    end

    respTime{condIdx} = respTimes;
    condIdx = condIdx+1;
end