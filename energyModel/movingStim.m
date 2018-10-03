clearvars;

%% create stimuli
% moving stim similar to the exp but 4 locations so that it is going to the
% right (net motion different from 0)
f1 = [ones(20,50) zeros(20,150)];
cycle4 = [f1; circshift(f1,50,2); circshift(f1,100,2); circshift(f1,150,2)];
stim = repmat(cycle4,2,1);
figure;imagesc(stim)

%% get the filters (general)
[left_1, left_2, right_1, right_2, dx, dt] = createFilters;

%% next steps depend on stim
% Step 3a: Define the space and time dimensions of the stimulus
% SPACE: x_stim is a row vector to hold sampled x-positions of the space.
stim_width=4;  %half width in degrees, gives 8 degrees total
x_stim=-stim_width:dx:round(stim_width-dx);

% TIME: t_stim is a col vector to hold sampled time intervals of the space
stim_dur=1.6;    %total duration of the stimulus in seconds
t_stim=(0:dt:round(stim_dur-dt))';


% Step 3c: convolve

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
%   Generate motion contrast matrix
energy_opponent = energy_right - energy_left; % L-R difference matrix
[xv yv] = size(energy_left); % Get the size of the response matrix
energy_flicker = total_energy/(xv * yv); % A value for average total energy

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

saveas(gcf,'realMotion.png','png' );

save('sumReal.mat','energy_flicker','motion_energy','total_energy')

