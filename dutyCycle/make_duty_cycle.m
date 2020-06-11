%% 
%Script to make figures for SSVEP review

%Change matlab default plotting options for better publication quality plots
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',30)
set(0,'DefaultTextFontSize',30)
myCol = [0 0 0; 0 0 0];
set(0,'DefaultAxesColorOrder',myCol)
set(0,'DefaultAxesLineStyleOrder',{'-',':'})
set(0,'DefaultAxesLineWidth',3)

%%
figure(400)
nPix = 256;
w=2;
[x y] = ndgrid(linspace(-pi,pi,nPix),linspace(-pi,pi,nPix));

stim = sin(w*y);

imagesc(stim)
colormap(gray)
axis off
%% Generate onset/offset response (i.e. odd and even harmonics)
nT=500;
t=linspace(0,1,nT);

% wList= [ 2 4 6 8 10];
% wAmp = complex([0 .707 0 1 0],[.5 .707 .5 0 .5]);     

 wList= [1:1:12]';
 nF = length(wList);
 wAmp = [linspace(.05,1,nF/2) linspace(1,.25,nF/2) ]';
 wAmp(2:2:end) = wAmp(2:2:end)*1;
 %wPhase = (-1*wList)+repmat([3*pi/2 4*pi/3]',nF/2,1);
 wPhase = (-1*wList)+repmat([0*pi/2 0*pi/2]',nF/2,1);
 
 wAmp = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;
 
%Make a fourier basis set for an easy waveform creator:
dutyCycle= .75
fourierBasis = dftmtx(nT);
invFourier  = conj(fourierBasis)/nT;
waveForm = wAmp'*fourierBasis(wList+1,:);
waveForm = real(waveForm);
waveForm = waveForm + circshift(waveForm,round(nT*dutyCycle));
waveForm = waveForm./max(abs(waveForm));
%

wAmp= 12*(abs(invFourier(wList+1,:)*waveForm'));


%timecourse panels
    
figure(300)
clf;

plot(t,waveForm,'k');
haxes1=gca;
hold on;
xlabel('Time (seconds)')
ylabel('Amplitude');
ylim([-1.2 1.2])
set(gca,'box','off')

% xticklocs = linspace(0,1,7);
% xticklabels = num2str(mod(xticklocs*360*w,360)','%3.0f');
% if i==1,
%     hfig1_pos = get(gcf,'Position')
% else
%     set(gcf,'Position',hfig1_pos);
% end

haxes1_pos = get(haxes1,'Position'); % store position of first axes
haxes1_pos(2) = haxes1_pos(2)*1;
haxes1_pos(2) = haxes1_pos(2)*1.1;
haxes1_pos(3) = haxes1_pos(3)*1;
haxes1_pos(4) = haxes1_pos(4)*.7;

set(haxes1,'Position',haxes1_pos);
% haxes2 = axes('Position',haxes1_pos,...
%     'XAxisLocation','top',...
%     'Color','none','YTick',[],'XTick',xticklocs,'XTickLabel',xticklabels);
% xlabel(haxes2,'Phase (degrees)')


%Make spec plot
figure(301);
clf;
%set(gcf,'Position',hfig1_pos);
freqs = 0:.5:20;
amps = 0*ones(size(freqs));
amps(2*wList+1)=abs(wAmp);

if exist('pdSpecPlot')
    pdSpecPlot(freqs,amps,[]);
else
    bar(freqs,amps);
end
% haxes1=gca;
% haxes1_pos = get(haxes1,'Position'); % store position of first axes
% haxes1_pos(2) = haxes1_pos(2)*1;
% haxes1_pos(2) = haxes1_pos(2)*1.1;
% haxes1_pos(3) = haxes1_pos(3)*1;
% haxes1_pos(4) = haxes1_pos(4)*.8;

%set(haxes1,'Position',haxes1_pos);
%ylim([0 1])
xlim([0 15])
xlabel('Frequency (Hertz)')
ylabel('Amplitude');
set(gca,'box','off')

% 
% %% Generate reversal response (i.e. even harmonics)
% nT=500;
% t=linspace(0,1,nT);
% 
% % wList= [ 2 4 6 8 10];
% % wAmp = complex([0 .707 0 1 0],[.5 .707 .5 0 .5]);     
%  
% wList= [1:1:12]';
%  nF = length(wList);
%  wAmp = [linspace(.05,1,nF/2) linspace(1,.25,nF/2) ]';
% wAmp(1:2:end) = 0;
%  wPhase = (-1*wList)+repmat([0 pi/2]',nF/2,1);
%  
%  wAmp = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;
%  
% %Make a fourier basis set for an easy waveform creator:
% 
% fourierBasis = dftmtx(nT);
% waveForm = wAmp'*fourierBasis(wList+1,:);
% waveForm = real(waveForm);
% waveForm = waveForm./max(abs(waveForm));
% %
% 
% %timecourse panels
%     
% figure(310)
% clf;
% 
% plot(t,waveForm,'k');
% haxes1=gca;
% hold on;
% xlabel('Time (seconds)')
% ylabel('Amplitude');
% ylim([-1.2 1.2])
% set(gca,'box','off')
% 
% xticklocs = linspace(0,1,7);
% xticklabels = num2str(mod(xticklocs*360*w,360)','%3.0f');
% if i==1,
%     hfig1_pos = get(gcf,'Position')
% else
%     set(gcf,'Position',hfig1_pos);
% 
% end
% haxes1_pos = get(haxes1,'Position'); % store position of first axes
% haxes1_pos(2) = haxes1_pos(2)*1;
% haxes1_pos(2) = haxes1_pos(2)*1.1;
% haxes1_pos(3) = haxes1_pos(3)*1;
% haxes1_pos(4) = haxes1_pos(4)*.7;
% 
% set(haxes1,'Position',haxes1_pos);
% % haxes2 = axes('Position',haxes1_pos,...
% %     'XAxisLocation','top',...
% %     'Color','none','YTick',[],'XTick',xticklocs,'XTickLabel',xticklabels);
% % xlabel(haxes2,'Phase (degrees)')
% 
% 
% %Make spec plot
% figure(311);
% clf;
% %set(gcf,'Position',hfig1_pos);
% freqs = 0:.5:20;
% amps = 0*ones(size(freqs));
% amps(2*wList+1)=abs(wAmp);
% pdSpecPlot(freqs,amps,[]);
% 
% % haxes1=gca;
% % haxes1_pos = get(haxes1,'Position'); % store position of first axes
% % haxes1_pos(2) = haxes1_pos(2)*1;
% % haxes1_pos(2) = haxes1_pos(2)*1.1;
% % haxes1_pos(3) = haxes1_pos(3)*1;
% % haxes1_pos(4) = haxes1_pos(4)*.8;
% 
% %set(haxes1,'Position',haxes1_pos);
% ylim([0 1])
% xlim([0 15])
% xlabel('Frequency (Hertz)')
% ylabel('Amplitude');
% set(gca,'box','off')










