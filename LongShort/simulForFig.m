
load('C:\Users\Marlene\Documents\git\dataLR\LRshortDC\Axx\Axx_LRshortDCV2_S01.mat')

% spatial interaction
elec = 23;
y1 = Axx(8).wave(elec,:); 
y2 =  Axx(9).wave(elec,:); 
simultO = Axx(10).wave(elec,:);  
simultP = y1+y2; 
spaInt = simultO - simultP; 
figure; hold on; plot(y1);plot(y2); ylim([-15 15]);saveas(gcf,'figures/simY1Y2','pdf')
figure; hold on; plot(simultO); plot(simultP); plot(spaInt);ylim([-15 15]);saveas(gcf,'figures/simSpa','pdf')

% temporal interaction
y1shift = circshift(y1,[0 length(y1)/2]);
y2shift = circshift(y2,[0 length(y2)/2]);
tempLeftO = Axx(11).wave(elec,:);
tempRightO = Axx(12).wave(elec,:);
tempLeftP = y1 + y1shift;
tempRightP = y2 + y2shift;
tempLeftInt = tempLeftO - tempLeftP;
tempRightInt = tempRightO - tempRightP;
figure; hold on;  plot(y1);  plot(y1shift); ylim([-15 15]); saveas(gcf,'figures/simY1Y1shift','pdf')
figure; hold on; plot(tempLeftO); plot(tempLeftP); plot(tempLeftInt); ylim([-15 15]);saveas(gcf,'figures/simTemp','pdf')


% linear prediction
figure; hold on;  plot(y1+y2shift);

% linear + spatial
spaIntshift = circshift(spaInt,[0 length(spaInt)/2]);
plot((y1+spaInt/2) + (y2shift+spaIntshift/2));

% linear + temporal
tempRightIntShift = circshift(tempRightInt,[0 length(tempRightInt)/2]);
plot((y1+tempLeftInt/2) + (y2shift+tempRightIntShift/2));

% linear + temporal + spatial
plot((y1+tempLeftInt/2+spaInt/2) + (y2shift+tempRightIntShift/2+spaIntshift/2));
ylim([-15 15]);
saveas(gcf,'figures/simPredictions','pdf')


% signal that is summed for the predictions
figure; hold on;  plot(spaInt/2); plot(spaIntshift/2);
plot(y1); plot(y2shift); plot(tempLeftInt/2); plot(tempRightIntShift/2);
ylim([-15 15]);
legend('S','Sshift','Y1','Y2shift','Tl','Tr')
saveas(gcf,'figures/simSig','pdf')