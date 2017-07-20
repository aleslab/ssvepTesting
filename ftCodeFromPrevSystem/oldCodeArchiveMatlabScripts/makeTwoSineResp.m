function nResp = makeTwoSineResp(c1, c2)

% takes two contrast values
% generate the response of the normalization model

% according to Carandini et al (1997)
% assumes linear receptive field (possibly with static nonlinearity)
% followed by divisive inhibition

% parameters of the normalization model are set here.

rMax = 1;
c50  = 13;    % contrast in percent rather than decimal
%contrast at which the response is half the maximum
orderNum = 2;
orderDen = 2;
base = 0;

% sinusoid parameters

w1 = 5;
w2 = 7;

ph2 = 0;

% Make sinusoids
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector


% raised cosine to simulate ON-OFF flickering stimuli
% input{1} = c1*(sin(2*pi*w1*t)) + c1;
% input{2} = c2*(sin(2*pi*w2*t+ph2)) + c2;

% sine to simulate contrast reversing
input{1} = c1*(sin(2*pi*w1*t)) ;
input{2} = c2*(sin(2*pi*w2*t+ph2)) ;

% input{1} = HalfRectify(input{1});
% input{2} = HalfRectify(input{2});

% input{1} = abs(input{1});
% input{2} = abs(input{2});


% Make output 

nResp = rMax*(input{1}+input{2}).^orderNum  ./ ...  % numerator
       ((input{1}+input{2}).^orderDen + c50.^orderDen) + base;  %denominator


% DEBUG code; comment out when not needed

figure(100)
clf
nRows = 3;

subplot(nRows,1,1)
plot(input{1},'k','linewidth',4)
hold on
plot(input{2},'--k','linewidth',4)
subplot(nRows, 2, 1)
plot(nResp, '--k', 'linewidth', 4)


title('Input','fontsize',20,'fontname','arial')
% axis off
grid on

% END DEBUG