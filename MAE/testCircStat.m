% test CircStat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orientationsAngle = [40 90 200; 10 7 100; 40 80 210; 10 10 150]; % ori in ex2
orientations = deg2rad(orientationsAngle);
amplitudes = [1 0.3 1.5; 0.8 0.5 1; 1 0.8 0.5; 0.2 0.2 1]; %w in ex2
amplitudes = amplitudes*10;
% 4 sbj, 3 conditions

for j = 1:3
  % confidence limits on mean angle
  confMean(j) = circ_confmean(orientations(:,j),[],amplitudes(:,j));
end

figure
for j = 1:3
  subplot(1,3,j);hold on;
  
  % compute and plot mean resultant vector length and direction
  r = circ_r(orientations(:,j));
  phi = circ_mean(orientations(:,j));
  zm = r*exp(i*phi');
  plot([0 real(zm)], [0, imag(zm)],'r','linewidth',1.5)
  
%   % plot the tuning function of the three neurons 
%   polar([orientations(:,j)' orientations(1,j)], [amplitudes(:,j)' amplitudes(1,j)],'k')
  
  % draw a unit circle
  zz = exp(i*linspace(0, 2*pi, 101)) * mw;
  plot(real(zz),imag(zz),'k:')
  plot([-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')

  axis square
  axis([-mw mw -mw mw])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  title(['condition ' num2str(j)])
end

figure
for j = 1:3
  subplot(1,3,j);hold on;
  
  % compute and plot mean resultant vector length and direction
  mw = max(max(amplitudes));
  r = circ_r(orientations(:,j),amplitudes(:,j)) * mw;
  phi = circ_mean(orientations(:,j),amplitudes(:,j));
  zm = r*exp(i*phi');
  plot([0 real(zm)], [0, imag(zm)],'r','linewidth',1.5)
  
%   % plot the tuning function of the three neurons 
%   polar([orientations(:,j)' orientations(1,j)], [amplitudes(:,j)' amplitudes(1,j)],'k')
  
  % draw a unit circle
  zz = exp(i*linspace(0, 2*pi, 101)) * mw;
  plot(real(zz),imag(zz),'k:')
  plot([-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')

  axis square
  axis([-mw mw -mw mw])
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  title(['condition ' num2str(j)])
end




  
color={'b','r','y'};
figure;hold on
mw = max(max(amplitudes));
% draw a unit circle
zz = exp(i*linspace(0, 2*pi, 101)) * mw;
plot(real(zz),imag(zz),'k:')
plot([-mw mw], [0 0], 'k:', [0 0], [-mw mw], 'k:')
axis square
% axis([-mw mw -mw mw])
set(gca,'xtick',[])
set(gca,'ytick',[])
for j = 1:3
    % compute and plot mean resultant vector length and direction
    r = circ_r(orientations(:,j),amplitudes(:,j)) * mw;
    phi = circ_mean(orientations(:,j),amplitudes(:,j));
    zm = r*exp(i*phi');
    plot([0 real(zm)], [0, imag(zm)],color{j},'linewidth',1.5)
    % should plot the variance here somehow... 
    polar([orientations(:,j)' orientations(1,j)], [amplitudes(:,j)' amplitudes(1,j)],'k')
end
color={'b','r','y'};
figure;hold on
zz = exp(i*linspace(0, 2*pi, 101));
plot(real(zz),imag(zz),'k:')
axis square
set(gca,'xtick',[])
set(gca,'ytick',[])
for j = 1:3
    % compute and plot mean resultant vector length and direction
    r = circ_r(orientations(:,j));
    phi = circ_mean(orientations(:,j));
    zm = r*exp(i*phi');
    plot([0 real(zm)], [0, imag(zm)],color{j},'linewidth',1.5)
end