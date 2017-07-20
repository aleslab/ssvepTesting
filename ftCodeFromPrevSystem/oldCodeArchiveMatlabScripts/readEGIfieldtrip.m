%%readEGIfieldtrip
clear all

%correct for powerdivaVideo sampling rate
timingCorrection = 432/1000;


%read in dataset
cfg = [];
cfg.dataset = '/Volumes/MRI/data/4D2/cptAmblyopes/fieldtrip/skeri0068/Raw/DF_amblyCPT4_12042009.raw';
elecfile = '/Volumes/MRI/data/4D2/cptAmblyopes/fieldtrip/skeri0068/Polhemus/DF_AmblyCPT4_12042009.elp';


[data] = ft_preprocessing(cfg);


cfg.trialdef.eventtype = '?';
cfg.trialdef.eventvalue = '?';

[cfg1] = ft_definetrial(cfg);

condList = {cfg1.event(:).value};
condList = unique(condList(2:end));
results = regexp(condList, 'c\d\d\d');
resIdx = ~cellfun('isempty', results);


condList = condList(resIdx);
clear results; clear resIdx; clear cfg1;





%%


cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = condList;
cfg.trialdef.prestim = 0 * timingCorrection;
cfg.trialdef.poststim = 2 * timingCorrection;


[cfg] = ft_definetrial(cfg);

elp = readelp(elecfile);

temppnt = [elp.X; elp.Y; elp.Z]';
elec.pnt(1:128,:) = temppnt(5:end,:);
elec.pnt(129,:) = temppnt(4,:);

templabel = {elp.labels};
elec.label = data.label;
% elec.label = templabel(5:end);
% elec.label{129} = templabel{4};

cfg.elec = elec;
clear elp; clear tmppnt; clear templabel; clear elec; clear elecfile; clear data;


%%




%artifact labeling
cfg.method = 'summary';
cfg.artfctdef.eog.channel = 'e17';
cfg.artfctdef.eog.feedback = 'yes';

[cfg, artifact] = ft_artifact_eog(cfg, data);

% cfg1.artfctdef.muscle.channel =  { 'e6' 'e7' 'e70' 'e71' 'e72' 'e73' 'e74' 'e75' 'e76' 'e77' 'e78' 'e79' 'e80' 'e81' 'e82' 'e83' 'e84' 'e85' 'e86' 'e87' };
% cfg1.artfctdef.muscle.feedback = 'yes';
%[cfg1, artifact] = ft_artifact_muscle(cfg1, data);


%%
%trial definition
cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = 'c003';
cfg.trialdef.prestim = 0 * timingCorrection;
cfg.trialdef.poststim = 2 * timingCorrection;
cfg.trialfun = 'lock2din';

[cfg] = ft_definetrial(cfg);
[trlData] = ft_redefinetrial(cfg, data);


cfg.artfctdef.reject = 'partial';
cleanData = ft_rejectartifact(cfg, trlData);


cfg.keepchannel = 'nan';
[cleanData] = ft_rejectvisual(cfg, cleanData);

numChans = length(cleanData.label);
badchannels = [];
for i = 1:numChans
	if isnan(cleanData.trial{1}(i,1))
		badchannels = [badchannels, cleanData.label(i,1)];
	end
end

cfg.badchannel = badchannels;
cleanData.elec = cfg.elec;
cleanData.feedback = 'yes';
cfg.feedback = 'yes';

[interp] = ft_channelrepair(cfg,cleanData);

cfg.reref = 'yes';
cfg.channel = 'all';
cfg.refchannel = 'all';
cfg.trackconfig = 'on';

[trlData] = ft_redefinetrial(cfg, trlData);
[rawAvgData] = ft_timelockanalysis(cfg, trlData);


[cleanData] = ft_redefinetrial(cfg, interp);






%filtering, segmenting and artifact rejection
cfg.lpfilter = 'yes';
cfg.lpfreq = 60;



layout = ft_prepare_layout(cfg, data);
layout.label = data.label;
cfg.layout = layout;
clear layout;



cfg.artfctdef.reject = 'partial';
%cleanData = ft_rejectartifact(cfg1, trlData);

%timelock averaging
cfg.keeptrials = 'yes';
cfg.vartrllength = 1;
[cleanAvgData] = ft_timelockanalysis(cfg, cleanData);

figure;
plot(rawAvgData.time/timingCorrection, rawAvgData.avg(:,:));
figure;
plot(cleanAvgData.time/timingCorrection, cleanAvgData.avg(:,:));
figure;
plot(rawAvgData.time/timingCorrection, rawAvgData.avg(6,:));
figure;
plot(cleanAvgData.time/timingCorrection, cleanAvgData.avg(6,:));






