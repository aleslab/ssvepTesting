function [cfg, cleanData, chanSubs, interpData] = ft_channelThreshold(cfg, data)
%Searches by trial for channels with amplitudes above a given threshold and
%then interpolates those channels for that trial and that trial alone

%50uv is arbitrary threshold for ERP data
chanThresh = 50;
numTrials =  size(data.trial, 2);
chanSubs = {};
cleanData = data;
for trlOp = 1: numTrials
	tmpcfg = cfg;
	tmpcfg.trials = trlOp;
    tmpcfg.badchannel = {};
    
	
	for nCh = 1:128
		tmpArray(nCh,trlOp) = max(abs(data.trial{trlOp}(nCh,:)));
		
		if  tmpArray(nCh,trlOp) > chanThresh;
			tmpcfg.badchannel = [tmpcfg.badchannel data.label{nCh}];
			
        end
        
	end
    
%if there are more than 8 (arbitrary amount) bad channels in a trial, Reset
%the threshold to 100 and detect again
    if numel(tmpcfg.badchannel) > 8
       tmpcfg.badchannel = {};
       for nCh = 1:128
           if tmpArray(nCh, trlOp) > 100;
               tmpcfg.badchannel = [tmpcfg.badchannel data.label{nCh}];
           end
       end
    end   
    chanSubs{trlOp} = tmpcfg.badchannel;
	
	
%interpolate bad channels    
    if  ~isempty(tmpcfg.badchannel)
        trlOp
        interpData = ft_channelrepair(tmpcfg, data);
        
        cleanData.trial{trlOp} = interpData.trial{1};
    end
    
 
end

%keeps track of bad channels
cleanData.chanSubs = chanSubs;


cfg = cleanData.cfg;

%corrects condition labels
cfg.eventold = cfg.event;
cfg.event = cell(1,size(cfg.trl,1));

newTrials = size(cfg.trl,1);

for n = 1:newTrials
	
	oldTrialNum = find((cfg.trl(n,1) >= cfg.trlold(:,1)) & (cfg.trl(n,1) <= cfg.trlold(:,2)));
	cfg.event{n} = cfg.eventold{oldTrialNum};
		

end


