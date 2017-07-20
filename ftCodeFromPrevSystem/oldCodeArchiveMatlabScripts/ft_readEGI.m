function [cfg, data, condList] = ft_readEGI(dataset, polhemus)



timingCorrection = 432/1000;
cfg = [];
%cfg.dataset = '/Volumes/MRI/data/4D2/cptAmblyopes/fieldtrip/skeri0068/Raw/DF_amblyCPT4_12042009.raw';
%elecfile = '/Volumes/MRI/data/4D2/cptAmblyopes/fieldtrip/skeri0068/Polhemus/DF_AmblyCPT4_12042009.elp';
cfg.dataset = dataset;
elecfile = polhemus;


cfg.trialdef.eventtype = '?';
cfg.trialdef.eventvalue = '?';

[cfg1] = ft_definetrial(cfg);

condList = {cfg1.event(:).value};
condList = unique(condList(2:end));
results = regexp(condList, 'c\d\d\d');
resIdx = ~cellfun('isempty', results);


condList = condList(resIdx);
clear results; clear resIdx; clear cfg1;



cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = condList;
cfg.trialdef.prestim = 0 * timingCorrection;
cfg.trialdef.poststim = 2 * timingCorrection;
cfg.trialfun = 'lock2dinCPT';


[cfg] = ft_definetrial(cfg);

elp  = readelp(elecfile);

temppnt = [elp.X; elp.Y; elp.Z]';
elec.pnt(1:128,:) = temppnt(5:end,:);
elec.pnt(129,:) = temppnt(4,:);

templabel = strcat('e', {elp.labels});

elec.label = templabel(5:end);
elec.label{129} = templabel{4};

cfg.elec = elec;
clear elp; clear tmppnt; clear templabel; clear elec; clear elecfile;

cfg.artfctdef.eog.trlpadding   = 0.1;
cfg.artfctdef.eog.fltpadding   = 0.1;
cfg.artfctdef.eog.artpadding   = 0.1;
cfg.artfctdef.eog.channel = 'e17';
cfg.artfctdef.eog.feedback = 'yes';

cfg.eventstart = cfg.trl(:,1,1);

[cfg, artifact] = ft_artifact_eog(cfg);

cfg.artfctdef.reject = 'partial';
[cfg] = ft_rejectartifact(cfg);

cfg.eventold = cfg.event;
cfg.event = cell(1,size(cfg.trl,1));

for nNew=1:size(cfg.trl,1)
	
	oldTrialNum = find((cfg.trl(nNew,1) >= cfg.trlold(:,1)) & (cfg.trl(nNew,1) <= cfg.trlold(:,2)));
	cfg.event{nNew} = cfg.eventold{oldTrialNum};
		

end
			
[data] = ft_preprocessing(cfg);

cfg.method = 'summary';
cfg.keepchannel = 'nan';
data.elec = cfg.elec;



[cleanData] = ft_rejectvisual(cfg, data);

numChans = length(cleanData.label);
badchannels = [];
for i = 1:numChans
	if isnan(cleanData.trial{1}(i,1))
		badchannels = [badchannels, cleanData.label(i,1)];
	end
end

cfg = cleanData.cfg;

cfg.eventold = cfg.event;
cfg.event = cell(1,size(cfg.trl,1));

newTrials = size(cfg.trl,1);

for n = 1:newTrials
	
	oldTrialNum = find((cfg.trl(n,1) >= cfg.trlold(:,1)) & (cfg.trl(n,1) <= cfg.trlold(:,2)));
	cfg.event{n} = cfg.eventold{oldTrialNum};
		

end

cfg = rmfield(cfg, 'method');

cleanData.feedback = 'yes';
cfg.feedback = 'yes';
cfg.reref = 'yes';
numChannels = cleanData.label;
numChannels(129)= [];
cfg.channel = numChannels;
cfg.refchannel = numChannels;
cfg.trackconfig = 'on';

[data] = ft_preprocessing(cfg);

cfg.badchannel = badchannels;
cfg.elec = cleanData.elec;
data.elec = cleanData.elec;
% 
% cfg.eventold = cfg.event;
% cfg.event = cell(1,size(cfg.trl,1));
% 
% newTrials = size(cfg.trl,1);
% 
% for n = 1:newTrials
% 	
% 	oldTrialNum = find((cfg.trl(n,1) >= cfg.trlold(:,1)) & (cfg.trl(n,1) <= cfg.trlold(:,2)));
% 	cfg.event{n} = cfg.eventold{oldTrialNum};
% 		
% 
% end


layout = ft_prepare_layout(cfg, data);
layout.label = data.label;
cfg.layout = layout;
clear layout; clear newTrials; 

% cfg.reref = 'yes';
% cfg.channel = 'all';
% cfg.refchannel = 'all';
% cfg.trackconfig = 'on';


[data] = ft_channelrepair(cfg,data);
%[data] = ft_redefinetrial(cfg, data);


cfg.keeptrials = 'yes';
cfg.vartrllength = 1;
data.cfg = cfg;


