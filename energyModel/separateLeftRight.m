clearvars
plotting = 0;

% filters parameters:
% width of temporal filters doesn't influence much
% beta does change energy for moving stim ( bell vs psychometric curve)
% width spatial filters only change scale


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERS PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1a: Define the space axis of the filters
nx=80;              %Number of spatial samples in the filter
max_x =2.0;         %2.0 % Half-width of filter (deg)
dx = (max_x*2)/nx;  %Spatial sampling interval of filter (deg)

% Spatial filter parameters
sx=0.5;   %0.5 % standard deviation of Gaussian, in deg.
sf=1.1;  % 1.1 % spatial frequency of carrier, in cpd

% Step 1b: Define the time axis of the filters
nt=100;         % Number of temporal samples in the filter
max_t=0.5;      % 0.5 Duration of impulse response (sec)
dt = max_t/nt;  % Temporal sampling interval (sec)

% Temporal filter parameters
k = 100;    % Scales the response into time units
slow_n = 9; % 9 % Width of the slow temporal filter
fast_n = 6; % 6 % Width of the fast temporal filter
beta = 0.9;  % 0.9
            % Beta. Represents the weighting of the negative
            % phase of the temporal relative to the positive 
            % phase.

            
% A row vector holding spatial sampling intervals
x_filt=linspace(-max_x,max_x,nx);
% A column vector holding temporal sampling intervals
t_filt=linspace(0,max_t,nt)';

% Spatial filter response
gauss=exp(-x_filt.^2/sx.^2);          %Gaussian envelope
even_x=cos(2*pi*sf*x_filt).*gauss;   %Even Gabor
odd_x=sin(2*pi*sf*x_filt).*gauss;    %Odd Gabor            
            
            
% Temporal filter response (formula as in Adelson & Bergen, 1985, Eq. 1)
slow_t=(k*t_filt).^slow_n .* exp(-k*t_filt).*(1/factorial(slow_n)-beta.*((k*t_filt).^2)/factorial(slow_n+2));
fast_t=(k*t_filt).^fast_n .* exp(-k*t_filt).*(1/factorial(fast_n)-beta.*((k*t_filt).^2)/factorial(fast_n+2));

% Step 1c: combine space and time to create spatiotemporal filters
e_slow= slow_t * even_x;    %SE/TS
e_fast= fast_t * even_x ;   %SE/TF
o_slow = slow_t * odd_x ;   %SO/TS
o_fast = fast_t * odd_x ;   % SO/TF

%--------------------------------------------------------------------------
%         STEP 2: Create spatiotemporally oriented filters
%--------------------------------------------------------------------------
left_1=o_fast+e_slow;      % L1
left_2=-o_slow+e_fast;     % L2
right_1=-o_fast+e_slow;    % R1
right_2=o_slow+e_fast;     % R2            
            
            
% figure;hold on;plot(slow_t);plot(fast_t)
% figure;imagesc(left_1)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STIMULUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3a: Define the space and time dimensions of the stimulus
% stimulus spatial resolution:
% 160 width for 8 degree
% 1 degree = 20 steps, 1 step = 0.05 deg

% stimulus temporal resolution:
% 200 length for 1.6 sec
% 1 step = 8ms, 24 steps = 192 ms

% SPACE: x_stim is a row vector to hold sampled x-positions of the space.
stim_width=4;  %half width in degrees, gives 8 degrees total
x_stim=-stim_width:dx:round(stim_width-dx);

% TIME: t_stim is a col vector to hold sampled time intervals of the space
stim_dur=4;    %total duration of the stimulus in seconds
t_stim=(0:dt:stim_dur-dt)';


cc=0;
for cycles = [10 5 2.5]
    cc=cc+1;dd=0;
for dc=5:5:95 % [10 50 90] % [10 25 50 75 90]
    clear stim stimL static
dd=dd+1;
% static stimulus
stim = ((square(cycles*pi*t_stim,dc)) +1) ./2;  % square(t,duty)
% 5pi = 5 cycles in 1 sec, 200 ms ON, 200 OFF 
stimL = repmat(stim,1,20) ;
static = [zeros(length(stimL),70), stimL, zeros(length(stimL),70)];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3c: convolve stim with filters

% Rightward responses
resp_right_1=conv2(static,right_1,'full');
resp_right_2=conv2(static,right_2,'full');

% Leftward responses
resp_left_1=conv2(static,left_1,'full');
resp_left_2=conv2(static,left_2,'full');

%--------------------------------------------------------------------------
%         STEP 4: Square the filter output
%--------------------------------------------------------------------------
resp_left_1 = resp_left_1.^2;
resp_left_2 = resp_left_2.^2;
resp_right_1 = resp_right_1.^2;
resp_right_2 = resp_right_2.^2;

% sum and diff responses
energy_right = resp_right_1 + resp_right_2;
energy_left = resp_left_1 + resp_left_2;


if dc==50
    figure;
    subplot(2,2,1);imagesc(static)
    subplot(2,2,3); imagesc(energy_right+energy_left); colorbar;
    title('Energy (Right+Left)')
    subplot(2,2,4); imagesc(energy_right-energy_left);colorbar
    title('Opponent (Right-Left)')
    
    % pool across space
    subplot(2,2,2);hold on;
    plot(sum(energy_right+energy_left,2))
    plot(sum(energy_right-energy_left,2))
    legend('sum','diff')
    saveas(gcf,['StaticDetails' num2str(cycles) '.png'])
end

% peak to peak amplitude excluding borders
minB = length(slow_t); maxB = length(stim)-length(slow_t);
maxEnergy(cc,dd,1) = max(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
avEnergy(cc,dd,1) = mean(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
amp(cc,dd,1) = max(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2)) ...
    - min(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
amp_mot(cc,dd,1) = max(sum(energy_right(minB:maxB,:)-energy_left(minB:maxB,:),2)) ...
    - min(sum(energy_right(minB:maxB,:)-energy_left(minB:maxB,:),2));
end
end

% figure; imagesc(amp); colorbar; title('Energy staticStim')
% figure; imagesc(amp_mot); colorbar; title('Motion staticStim')



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOVING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeStim = 20; cc=0;
for cycles = [10 5 2.5]
    cc=cc+1;dd=0;
for dc=5:5:95 % [10 50 90]
    clear stim stimL stimR
    dd=dd+1;
stimL = ((square(cycles/2*pi*t_stim,dc/2)) +1) ./2;  % /2 for motion

stimR(1:length(stimL)/(cycles*2)) = 0;
stimR(length(stimL)/(cycles*2)+1:length(stimL)) = stimL(1:length(stimL)-length(stimL)/(cycles*2));
stimL = repmat(stimL,1,20) ;
stimR = repmat(stimR,20,1) ;
mouv = [zeros(length(stimL),57), stimL, zeros(length(stimL),6), stimR', zeros(length(stimR),57)];

% figure; hold on; plot(stimR');plot(stimL)
% plot(((square(cycles*pi*t_stim,dc)) +1) ./2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3c: convolve stim with filters

% Rightward responses
resp_right_1=conv2(mouv,right_1,'full');
resp_right_2=conv2(mouv,right_2,'full');

% Leftward responses
resp_left_1=conv2(mouv,left_1,'full');
resp_left_2=conv2(mouv,left_2,'full');

%--------------------------------------------------------------------------
%         STEP 4: Square the filter output
%--------------------------------------------------------------------------
resp_left_1 = resp_left_1.^2;
resp_left_2 = resp_left_2.^2;
resp_right_1 = resp_right_1.^2;
resp_right_2 = resp_right_2.^2;

% sum and diff responses
energy_right = resp_right_1 + resp_right_2;
energy_left = resp_left_1 + resp_left_2;

% figure
if dc == 50
figure;
subplot(2,2,1);imagesc(mouv)
subplot(2,2,3); imagesc(energy_right+energy_left); colorbar;
title('Energy (Right+Left)')
subplot(2,2,4); imagesc(energy_right-energy_left);colorbar
title('Opponent (Right-Left)')

% pool across space
subplot(2,2,2);hold on;
plot(sum(energy_right+energy_left,2))
plot(sum(energy_right-energy_left,2))
legend('sum','diff')
    saveas(gcf,['MotionDetails' num2str(cycles) '.png'])
end


% peak to peak amplitude excluding borders
minB = length(slow_t); maxB = length(mouv)-length(slow_t);
amp(cc,dd,2) = max(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2)) ...
    - min(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
amp_mot(cc,dd,2) = max(sum(energy_right(minB:maxB,:)-energy_left(minB:maxB,:),2)) ...
    - min(sum(energy_right(minB:maxB,:)-energy_left(minB:maxB,:),2));
maxEnergy(cc,dd,2) = max(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
avEnergy(cc,dd,2) = mean(sum(energy_right(minB:maxB,:)+energy_left(minB:maxB,:),2));
end
end

% figure; imagesc(ampM); colorbar;title('Energy movingStim')
% figure; imagesc(amp_mot_M); colorbar;title('motion movingStim')


%% PLOT
for mm =1:2
figure;hold on;
subplot(2,2,1);hold on;
for cc=1:3
    plot(amp(cc,:,mm),'LineWidth',2)
end
legend('10','5','2.5','Location','best')
title('Energy Amplitude')
subplot(2,2,2);hold on;
for cc=1:3
    plot(amp_mot(cc,:,mm),'LineWidth',2)
end
legend('10','5','2.5','Location','best')
title('Motion Energy')

subplot(2,2,3);hold on;
for cc=1:3
    plot(maxEnergy(cc,:,mm),'LineWidth',2)
end
title('Max Energy')

subplot(2,2,4);hold on;
for cc=1:3
    plot(avEnergy(cc,:,mm),'LineWidth',2)
end
title('Av Energy')

if mm==1
    saveas(gcf,'Static.png')
else
    saveas(gcf,'Mouving.png')
end
% saveas(gcf,['Swidth' num2str(max_x) '-Sfreq' num2str(sf) '-Timp' num2str(max_t) '-Tbeta' num2str(beta) '-Twidth' num2str(slow_n) '.png'])
end
%%
