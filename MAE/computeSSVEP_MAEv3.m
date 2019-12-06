 

dataOut = '/Users/marleneponcet/Documents/data/MAEv3/Axx/';

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


% get the condition labels
phase=[10 10 90];
overlap=[0 1 1];
direction = {'none', 'left','right'};
condNb = 1;
for dir = 1:3
    for cc=1:3
        Axx(condNb).condLabel = ['HighSF' direction{dir} 'G' num2str(overlap(cc)+1) 'phase' num2str(phase(cc))] ;
        condNb = condNb+1;
    end
end
save([dataOut 'Axx_S0' num2str(ff)],'Axx','cfg')

% get the condition labels v3b
phase=[90 10 90];
overlap=[0 1 1];
direction = {'none', 'left','right'};
condNb = 1;
for dir = 1:3
    for cc=1:3
        Axx(condNb).condLabel = ['LowSF' direction{dir} 'G' num2str(overlap(cc)+1) 'phase' num2str(phase(cc))] ;
        condNb = condNb+1;
    end
end
save([dataOut 'Axx_S54'],'Axx','cfg')



cc = 1;
for aa=1:3
    [Axx(9+cc)] = substractAxx(Axx(aa+3),Axx(aa));
    Axx(9+cc).condLabel = ['sub-' Axx(aa+3).condLabel];
    [Axx(10+cc)] = substractAxx(Axx(aa+6),Axx(aa));
    Axx(10+cc).condLabel = ['sub-' Axx(aa+6).condLabel];    
    cc = cc+2;
end

save([dataOut 'Axx_S0' num2str(ff) 'sub'],'Axx','cfg')
