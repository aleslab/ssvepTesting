load('/Users/marleneponcet/Documents/dataLongRangeV2/AxxFiles/Axx_S05.mat')

signal = Axx.cond(1).wave(19,:); % get one channel
figure;plot(Axx.cond(1).time/1000,signal,'r');title('original signal');

nfr = floor(length(Axx.cond(1).time)/2)+1;

fftsig = fft(signal);
half = fftsig(1:nfr);
amplitude = abs(half);
dHz = 1000/(Axx.cond(1).time(end)+Axx.cond(1).dtms);
freq = 0: dHz : dHz*(nfr-1);
figure; bar(freq,amplitude,'r'); title('original signal');   

shift = circshift(signal,[0 length(signal)/2]); 
Motion.wave = signal + shift;
figure; plot(Axx.cond(1).time,shift,'g'); hold on; plot(Axx.cond(1).time,signal,'m'); plot(Axx.cond(1).time,Motion.wave,'b'); legend('shift','original','sum');
test = fft(Motion.wave);
half2 = test(2:floor(end/2));
figure;bar(freq,abs(half2));title('fft results')

periods = 1 ./(Axx.cond(1).freq);
phaseShift = Axx.cond(1).time(end/2+1)/1000 ./ periods .* 360; 
sinH = Axx.cond(1).sin(19,:);
cosH = Axx.cond(1).cos(19,:);
shiftSin = sind(phaseShift) .*  cosH + cosd(phaseShift) .* sinH;
shiftCos = cosd(phaseShift) .* cosH - sind(phaseShift) .* sinH ;

Motion.sin = sinH + shiftSin;
Motion.cos = cosH + shiftCos;
Motion.amp = sqrt(Motion.sin.^2 + Motion.cos.^2);
figure;bar(Motion.amp(1:100));title('cos-sin results')

figure; bar(Axx.cond(1).freq,Axx.cond(1).amp(19,:),'r');   

periods = 1 ./freq;
phaseShift = Axx.cond(1).time(end/2+1) ./ periods .* 360; 
sinH = Axx.cond(1).sin(19,:);
cosH = Axx.cond(1).cos(19,:);
shiftSin = sind(phaseShift) .*  cosH + cosd(phaseShift) .* sinH;
shiftCos = cosd(phaseShift) .* cosH - sind(phaseShift) .* sinH ;

Motion.sin = sinH + shiftSin;
Motion.cos = cosH + shiftCos;
Motion.amp = sqrt(Motion.sin.^2 + Motion.cos.^2);
figure;bar(Motion.amp(1:100));title('cos-sin results')
