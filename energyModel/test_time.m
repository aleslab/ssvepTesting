


% Step 1b: Define the time axis of the filters
nt=100;         % Number of temporal samples in the filter
max_t=0.5;      % Duration of impulse response (sec)
dt = max_t/nt;  % Temporal sampling interval (sec)

% A column vector holding temporal sampling intervals
t_filt=linspace(0,max_t,nt)';



% Temporal filter parameters
k = 80;    % Scales the response into time units
slow_n = 18; % Width of the slow temporal filter
fast_n = 12; % Width of the fast temporal filter
beta =1;  % Beta. Represents the weighting of the negative
            % phase of the temporal relative to the positive 
            % phase.

% Temporal filter response (formula as in Adelson & Bergen, 1985, Eq. 1)
slow_t=(k*t_filt).^slow_n .* exp(-k*t_filt).*(1/factorial(slow_n)-beta.*((k*t_filt).^2)/factorial(slow_n+2));
fast_t=(k*t_filt).^fast_n .* exp(-k*t_filt).*(1/factorial(fast_n)-beta.*((k*t_filt).^2)/factorial(fast_n+2));

figure;plot(slow_t+fast_t);hold on;plot(slow_t-fast_t)

t=0:dt:4;
stim = ((square(2*2*pi*t,20)) +1) ./2;


slow_out = conv(stim,slow_t,'full');
fast_out = conv(stim,fast_t,'full');

% figure; hold on; plot(slow_out);plot(fast_out);
% plot(stim)
% legend('slow','fast')

resp_slow = slow_out;
resp_fast = fast_out;

figure;
energy_right = resp_slow + resp_fast;
hold on;plot(energy_right,'LineWidth',2)
% plot(resp_slow,'LineWidth',2)
% plot(resp_fast,'LineWidth',2)
plot(resp_slow - resp_fast,'LineWidth',2)

