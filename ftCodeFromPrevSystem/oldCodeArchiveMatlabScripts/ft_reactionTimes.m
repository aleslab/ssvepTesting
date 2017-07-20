function [cfg] = ft_reactionTimes(projectInfo,optns)


inputDir = fullfile(projectInfo.projectDir,projectInfo.subjId,'Raw');
outputDir = projectInfo.powerDivaExportDir;


allRaw = dir(fullfile(inputDir,'*.raw'));
Raw = allRaw(length(allRaw)).name;


dataset = fullfile(inputDir, Raw);
elps = dir(fullfile(projectInfo.projectDir, projectInfo.subjId, 'Polhemus', '*.elp'));
polhemus = fullfile(projectInfo.projectDir, projectInfo.subjId, 'Polhemus', elps.name);

cfg =[];
cfg.dataset = dataset;
cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue = 'c003';

[trl, eventName, reactionTime] = responseLockCPT(cfg);


save(fullfile(projectInfo.projectDir,projectInfo.subjId, '_dev_','reactionTime.mat'),'eventName', 'reactionTime');