function [cfg, data, cond] = ft_responseLockEGI(dataset, polhemus) 

timingCorrection = 432/1000;
cfg = [];

% cfg.dataset = '/Volumes/MRI/data/4D2/CPT/test/skeri0087/Raw/SM_amblyCPT4_11182009.raw';
% elecfile = '/Volumes/MRI/data/4D2/CPT/test/skeri0087/Polhemus/SM_AmblyCPT4_11182009.elp'; 
cfg.dataset = dataset;
elecfile = polhemus;
cond = 'c007';

cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = cond;
cfg.trialdef.prestim = 1.5 * timingCorrection;
cfg.trialdef.poststim = .5 * timingCorrection;
cfg.trialfun = 'responseLockCPT';


[cfg] = ft_definetrial(cfg);

elp  = readelp(elecfile);
elec = [];
temppnt = [elp.X; elp.Y; elp.Z]';
elec.pnt(1:128,:) = temppnt(5:end,:);
elec.pnt(129,:) = temppnt(4,:);

templabel = strcat('e', {elp.labels});

elec.label = templabel(5:end);
elec.label{129} = templabel{4};

cfg.elec = elec;
clear elp; clear temppnt; clear templabel; clear elec; clear elecfile;


cfg.feedback = 'yes';
cfg.reref = 'yes';
numChannels = cfg.elec.label';
numChannels(129)= [];
cfg.channel = numChannels;
cfg.refchannel = numChannels;
cfg.trackconfig = 'on';


[data] = ft_preprocessing(cfg);
data.elec = cfg.elec;

cfg.method = 'summary';
cfg.keepchannel = 'nan';


[cleanData] = ft_rejectvisual(cfg, data);

numChans = length(cleanData.label);
badchannels = [];
for i = 1:numChans
	if isnan(cleanData.trial{1}(i,1))
		badchannels = [badchannels, cleanData.label(i,1)];
	end
end



cfg = cleanData.cfg;
cfg.badchannel = badchannels;
cfg.elec = data.elec;

[cleanData] = ft_channelrepair(cfg, cleanData);

cfg.artfctdef.eog.trlpadding   = 0.1;
cfg.artfctdef.eog.fltpadding   = 0.1;
cfg.artfctdef.eog.artpadding   = 0.1;
cfg.artfctdef.eog.channel = 'e17';
cfg.artfctdef.eog.feedback = 'yes';


[cfg, eogArtifact] = ft_artifact_eog(cfg);



cfg.artfctdef.reject = 'partial';
[cleanData] = ft_rejectartifact(cfg, cleanData);

[cfg, cleanData, chanSubs, interpData] = ft_channelThreshold(cfg, cleanData);



cfg = rmfield(cfg, 'method');




layout = ft_prepare_layout(cfg, cleanData);
layout.label = cleanData.label;
cfg.layout = layout;
clear layout; clear newTrials; 

%[data] = ft_redefinetrial(cfg, data);


cfg.keeptrials = 'yes';
cfg.vartrllength = 1;
cleanData.cfg = cfg;

data = cleanData;




% 
% cleanAvgData = ft_timelockanalysis(cfg, cleanData);
% rawAvgData = ft_timelockanalysis(cfg, data);
% 
% 
% figure
% plot(cleanAvgData.time/timingCorrection, cleanAvgData.avg(:,:));
% figure
% plot(rawAvgData.time/timingCorrection, rawAvgData.avg(:,:));
% figure
% plot(cleanAvgData.time/timingCorrection, cleanAvgData.avg(6,:));
% figure
% plot(rawAvgData.time/timingCorrection, rawAvgData.avg(6,:));