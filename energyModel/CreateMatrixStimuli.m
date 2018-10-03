

% create space-time matrix corresponding to the stimuli used in the Duty
% Cycle experiment
stimWidth = 8;

% try a moving stim with same dim as the real stim
% timeLength=1504ms = 2 left-right cycles in motion
small = diag(ones(1,stimWidth));
stim=[];
while length(stim)<100
    stim = [stim;small];
end
figure;imagesc(stim)
save('test','stim')


% STATIC conditions
totalCycle = [8/85 16/85 32/85]; % in Hz this is the onset of the single stimulus
onTime = [1/8 2/8 4/8 6/8 7/8];
% timeLength = round(32/85*1000);
minTotDuration = 6*32/85; % same as in the rating experiment
condition = 0;
for tc=1:length(totalCycle)
    for ot = 1:length(onTime)
        fullcycle = [];
        condition = condition+1;
        stimON = round(onTime(ot)*totalCycle(tc)*1000);
        stimOFF = round(totalCycle(tc)*1000 - stimON);
        cycle = [repmat(1,stimON,1); repmat(0,stimOFF,1)];
        while length(fullcycle)<minTotDuration
            fullcycle = [fullcycle;cycle];
        end
        stim = zeros(length(fullcycle),stimWidth);
        stim(:,3) = fullcycle;
        stim(:,4) = fullcycle;
        figure;imagesc(stim)
        save(['static' num2str(condition,'%.2d') '.mat'],'stim');
    end
end


% MOTION conditions
totalCycle = [8/85 16/85 32/85]; 
onTime = [1/8 2/8 4/8 6/8 7/8];
% timeLength = round(32/85*1000);
minTotDuration = 6*32/85;
condition = 0;
for tc=1:length(totalCycle)
    for ot = 1:length(onTime)
        fullcycleR = []; fullcycleL = [];
        condition = condition+1;
        stimON = round(onTime(ot)*totalCycle(tc)*1000);
        stimOFF = round(totalCycle(tc)*1000 - stimON);
        cycleRight = [repmat(1,stimON,1); repmat(0,stimOFF,1); repmat(0,round(totalCycle(tc)*1000),1)];
        cycleLeft = [repmat(0,round(totalCycle(tc)*1000),1);repmat(1,stimON,1); repmat(0,stimOFF,1)];
        while length(fullcycleL)<minTotDuration
            fullcycleL = [fullcycleL;cycleLeft];
            fullcycleR = [fullcycleR;cycleRight];
        end
        stim = zeros(length(fullcycleL),stimWidth);
        stim(:,3) = fullcycleL;
        stim(:,4) = fullcycleL;
        stim(:,6) = fullcycleR;
        stim(:,7) = fullcycleR;
        figure;imagesc(stim)
%         save(['motion' num2str(condition,'%.2d') '.mat'],'stim');
    end
end
% conditions 1,2,4,5,6,10,12,14 were not presented in the experiment