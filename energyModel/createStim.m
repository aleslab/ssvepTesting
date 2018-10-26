clearvars;
% add other durations: 100, 200, 400 for around 2.5,5,10Hz
% this should be changed in the timeON loop (+ repmat, etc.)

%% create stimuli
% should be 200 x 160 
% 200 = 8 degrees, separation of less than a degree = 8+8 = 16 = 0.64
% Will have 2 stim each of 46 = 1.84 degrees
% 160 = 1.6 s, each step = 10 ms, half cycle (one side) = 20, full cycle (2
% sides = motion) = 40

stimWidth = 46;
totWidth = 200;
midSep = 8;

% location short-range
locS1(1,:) = stimWidth+1:totWidth/2-midSep; locS2(1,:) = totWidth/2+midSep+1:totWidth-stimWidth;
% location long-range
locS1(2,:) = 1:stimWidth; locS2(2,:) = totWidth-stimWidth+1:totWidth;
for tt=1:2 % 1=SR, 2=LR
    for timeON = 1:20
    oneCycle = zeros(40,200);
    oneCycle(1:timeON,locS1(tt,:)) = ones(timeON,stimWidth);
    oneCycle(21:20+timeON,locS2(tt,:)) = ones(timeON,stimWidth);
    allStim(:,:,timeON,tt) = repmat(oneCycle,4,1);
%     figure;imagesc(allStim(:,:,timeON,tt))
    end
end


%% get the filters (general)
[left_1, left_2, right_1, right_2, dx, dt] = energyFilters;

%% next steps depend on stim
% Step 3a: Define the space and time dimensions of the stimulus
% SPACE: x_stim is a row vector to hold sampled x-positions of the space.
stim_width=4;  %half width in degrees, gives 8 degrees total
x_stim=-stim_width:dx:round(stim_width-dx);

% TIME: t_stim is a col vector to hold sampled time intervals of the space
stim_dur=1.6;    %total duration of the stimulus in seconds
t_stim=(0:dt:round(stim_dur-dt))';


% Step 3c: convolve

%% add loop here for all the stim to test
for tt = 1:2
for ss=1:size(allStim,3)
stim = allStim(:,:,ss,tt);
% Rightward responses
resp_right_1=conv2(stim,right_1,'valid');
resp_right_2=conv2(stim,right_2,'valid');

% Leftward responses
resp_left_1=conv2(stim,left_1,'valid');
resp_left_2=conv2(stim,left_2,'valid');
%--------------------------------------------------------------------------
%         STEP 4: Square the filter output
%--------------------------------------------------------------------------
resp_left_1 = resp_left_1.^2;
resp_left_2 = resp_left_2.^2;
resp_right_1 = resp_right_1.^2;
resp_right_2 = resp_right_2.^2;
%--------------------------------------------------------------------------
%         STEP 5: Normalise the filter output
%--------------------------------------------------------------------------
% Calc left and right energy
energy_right= resp_right_1 + resp_right_2;
energy_left= resp_left_1 + resp_left_2;

%%%%% Calc total energy
total_energy = sum(sum(energy_right))+sum(sum(energy_left));
totEnergy(tt,ss) = total_energy;

% Normalise each directional o/p by total output
RR1 = sum(sum(resp_right_1))/total_energy;
RR2 = sum(sum(resp_right_2))/total_energy;
LR1 = sum(sum(resp_left_1))/total_energy;
LR2 = sum(sum(resp_left_2))/total_energy;

%         STEP 6: Sum the paired filters in each direction
right_Total = RR1+RR2;
left_Total = LR1+LR2;
%         STEP 7: Calculate net energy as the R-L difference
motion_energy = right_Total - left_Total;
motionEnergy(tt,ss) = motion_energy;

%   Generate motion contrast matrix
energy_opponent = energy_right - energy_left; % L-R difference matrix
[xv yv] = size(energy_left); % Get the size of the response matrix
energy_flicker = total_energy/(xv * yv); % A value for average total energy
flickerEnergy(tt,ss) = energy_flicker;

if ss == 1 || ss==size(allStim,3)/2 || ss==size(allStim,3)
% Plot the stimulus
figure; 
subplot(2,1,1);
imagesc(stim); 
colormap(gray);
axis off
caxis([0 1.0]);
axis equal
title('Stimulus');

% Plot the output:
% Re-scale (normalize) each pixel in the L-R matrix using average energy.
motion_contrast = energy_opponent/energy_flicker;

% Plot, scaling by max L or R value
mc_max = max(max(motion_contrast));
mc_min = min(min(motion_contrast));
if (abs(mc_max) > abs(mc_min))
    peak = abs(mc_max);
else
    peak = abs(mc_min);
end

subplot(2,1,2);
imagesc(motion_contrast); 
colormap(gray);
axis off
caxis([-peak peak]); colorbar;
axis equal
title('Normalised Motion Energy');

saveas(gcf,['stim' num2str(tt) 'dc' num2str(ss/size(allStim,3)*100) 'png'],'png' )
end


end
end


load('sumReal.mat')

figure; 
subplot(1,3,1);hold on;
plot(totEnergy(1,:))
plot(totEnergy(2,:))
plot(20,total_energy,'Marker','*')
legend('SR','LR','motion')
title('totEnergy')

subplot(1,3,2);hold on;
plot(flickerEnergy(1,:))
plot(flickerEnergy(2,:))
plot(20,energy_flicker,'Marker','*')
legend('SR','LR','motion')
title('flickerEnergy')

subplot(1,3,3);hold on;
plot(motionEnergy(1,:))
plot(motionEnergy(2,:))
plot(20,motion_energy,'Marker','*')
legend('SR','LR','motion')
title('motionEnergy')

saveas(gcf,'Energy summary','png' )
