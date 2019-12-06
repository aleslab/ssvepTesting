 
dataOut = '/Users/marleneponcet/Documents/data/MAE/cleanComp/';

% compute steady state
cleanData=rmfield(cleanData,'cfg'); % clear cfg so that the data does not become too heavy

% 1 Axx per sbj
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
cfg.layout = 'biosemi128.lay';
allcond = unique(cleanData.trialinfo(:,1));

for cond=1:length(allcond)
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond));
    [Axx(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
end

%%% separate first and second half of the test period
for cond=1:length(allcond)
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond) & cleanData.trialinfo(:,4)<4);
    [AxxFirst(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
end
for cond=1:length(allcond)
    cfg.trials = find(cleanData.trialinfo(:,1) == allcond(cond) & cleanData.trialinfo(:,4)>3);
    [AxxSecond(cond)] = ft_steadystateanalysis(cfg, cleanData); % the field elec is present only if a channel has been replaced
end

% get the condition labels
spatialFq = [2 0.5 2]; % standard fq 0.5 /4 or *4
testFq = round([85/10 85/20 85/10]);
fovea = [1 1 0];
direction = {'none', 'left','right'};
condNb = 1;
for dir = 1:3
    for cc=1:3
        Axx(condNb).condLabel = [direction{dir} num2str(spatialFq(cc)) 'fq' num2str(testFq(cc)) 'fov' num2str(fovea(cc))] ;
        condNb = condNb+1;
    end
end

save([dataOut 'Axx_comp_S01'],'Axx','cfg')



% get the condition labels (S02 S05)
spatialFq = [2 0.5 2]; % standard fq 0.5 /4 or *4
testFq = round([85/20 85/10 85/10]);
direction = {'none', 'left','right'};
condNb = 1;
for dir = 1:3
    for cc=1:3
        Axx(condNb).condLabel = [direction{dir} num2str(spatialFq(cc)) 'fq' num2str(testFq(cc))] ;
        condNb = condNb+1;
    end
end
save([dataOut 'Axx_comp_S05'],'Axx','cfg')
