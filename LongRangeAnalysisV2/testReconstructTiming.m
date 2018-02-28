%%
%sampling
t = 0 : 0.01 : 2-0.01;
%Frequency 
f1 = 10; f2=5; 
%Signal
signal = sin(f1*2*pi*t) + sin(f2*2*pi*t) + randn(size(t));
figure;plot(t,signal);

fftsig = fft(signal);
half = fftsig(2:floor(end/2));
amplitude = abs(half);
cosH = real(half);
sinH= -imag(half);
    
freqs = [1:length(amplitude)] * 1/(t(end)+t(2));
figure; bar(freqs,amplitude,'r'); 

shift = circshift(signal,[0 20]);
periods = 1./freqs;
phaseShift =  0.2 ./ periods .* 360;

shiftCos = cosd(phaseShift) .* cosH - sind(phaseShift) .* sinH;
shiftSin = sind(phaseShift) .* cosH + cosd(phaseShift) .* sinH;

Motion.wave = signal +shift;
Motion.sin = sinH + shiftSin;
Motion.cos = cosH + shiftCos;
Motion.amp = sqrt(Motion.sin.^2 + Motion.cos.^2);

figure;plot(Motion.wave); hold on; plot(signal,'r');
figure;bar(freqs,Motion.amp)

test = fft(Motion.wave);
half2 = test(2:floor(end/2));
figure;bar(freqs,abs(half2))
hold on;
bar(freqs,amplitude,0.3,'FaceColor','r')

% %%
% t = Axx.cond.freq;
% fftsig = fft(signal);
% half = fftsig(2:floor(end/2));
% amplitude = abs(half);
% cos = real(half);
% sin= -imag(half);
% freqs = [1:length(amplitude)] * 1/(t(end)-t(1));
% figure; bar(freqs,amplitude,'r'); 
% 
