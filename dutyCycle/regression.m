%%%% compute slopes, 95% CI, etc. for comparing conditions 

clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults
addpath polyparci/

% load individual data
% dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
dataDir = 'C:\Users\Marlene\Documents\JUSTIN\data\dutyCycle\Axx\';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7:11 13:15 17:18 20]; 
dcVal = [12.5 25 50 75 87.5];

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    dataSbj(:,ss) = Axx;
end

pickElec = [23];
elec = 1;
% there must be a better way...
f1Amp = zeros(length(dataSbj),length(keepSbj));
for ss=1:length(keepSbj)
    for cond=1:15 % static = 1st harmonic
        f1Amp(cond,ss) = dataSbj(cond,ss).amp(pickElec(elec),dataSbj(cond,ss).i1f1);
    end
    for cond=16:22 % moving = 2nd harmonic
        f1Amp(cond,ss) = dataSbj(cond,ss).amp(pickElec(elec),dataSbj(cond,ss).i1f1*2-1);
    end    
end
% just for checking: plot the ind data on which the regression will be done
figure;hold on;
subplot(2,2,1);hold on;
for cond = 1:5
    scatter(repmat(cond,length(keepSbj),1),f1Amp(cond,:)')
end
subplot(2,2,2);hold on;
for cond = 6:10
    scatter(repmat(cond,length(keepSbj),1),f1Amp(cond,:)')
end
subplot(2,2,3);hold on;
for cond = 11:15
    scatter(repmat(cond,length(keepSbj),1),f1Amp(cond,:)')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Effect of DC on F1 amplitude?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% x = repmat(dcVal,length(keepSbj),1);
% y = f1Amp(1:5,:)';
% mdl = fitlm(x(:),y(:));
% model = table2array(mdl.Coefficients); 
% slope(1) = model(2,1);
% % confInt = conflevel * standard error
% ci95(1,:) = tinv([0.025 0.975], mdl.NumObservations-1) * model(2,2);
% ci99(1,:) = tinv([0.005 0.995], mdl.NumObservations-1) * model(2,2);
% 
% y = f1Amp(6:10,:)';
% mdl = fitlm(x(:),y(:));
% model = table2array(mdl.Coefficients); 
% slope(2) = model(2,1);
% ci95(2,:) = tinv([0.025 0.975], mdl.NumObservations-1) * model(2,2);
% ci99(2,:) = tinv([0.005 0.995], mdl.NumObservations-1) * model(2,2);
% y = f1Amp(11:15,:)';
% mdl = fitlm(x(:),y(:));
% model = table2array(mdl.Coefficients); 
% slope(3) = model(2,1);
% ci95(3,:) = tinv([0.025 0.975], mdl.NumObservations-1) * model(2,2);
% ci99(3,:) = tinv([0.005 0.995], mdl.NumObservations-1) * model(2,2);

% x = repmat(dcVal,length(keepSbj),1);
% y = f1Amp(1:5,:)';
% [coef,errorEst]=polyfit(x(:),y(:),1);
% % delta is the standard error estimate
% [p,S,mu] = polyfit(x,y,n);  
% figure;
% subplot(2,2,1); hold on
% plot(x(:),y(:),'bo')
% plot(x(:),y_fit,'r-')
% plot(x(:),y_fit+2*delta,'m--',x(:),y_fit-2*delta,'m--')
% 95% prediction interval = y±2?
% 
% subplot(2,2,2); hold on
% y = f1Amp(6:10,:)';
% [coef,errorEst]=polyfit(x(:),y(:),1);
% [y_fit,delta] = polyval(coef,x(:),errorEst);
% plot(x(:),y(:),'bo')
% plot(x(:),y_fit,'r-')
% plot(x(:),y_fit+2*delta,'m--',x(:),y_fit-2*delta,'m--')
% 
% subplot(2,2,3); hold on
% y = f1Amp(11:15,:)';
% [coef,errorEst]=polyfit(x(:),y(:),1);
% [y_fit,delta] = polyval(coef,x(:),errorEst);
% plot(x(:),y(:),'bo')
% plot(x(:),y_fit,'r-')
% plot(x(:),y_fit+2*delta,'m--',x(:),y_fit-2*delta,'m--')


x = repmat(dcVal,length(keepSbj),1);
y = f1Amp(1:5,:)'; % 10hz static
data(1,:) = y(:);
y = f1Amp(6:10,:)'; % 5Hz static
data(2,:) = y(:);
y = f1Amp(11:15,:)'; % 2.5Hz static
data(3,:) = y(:);

for cond = 1: 3
    [coefV,errorEst]=polyfit(x(:)',data(cond,:),1);
    coef(:,cond)=coefV;
    tmp = polyparci(coefV,errorEst,0.95); % gives lower and upper CI
    ci95(:,cond) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)]; % only keep the diff
    tmp = polyparci(coefV,errorEst,0.99); 
    ci99(:,cond) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)];
end

% don't include condition with no SSVEP
yRej = f1Amp(2:5,:)'; % 10Hz static
xRej = repmat(dcVal(2:5),length(keepSbj),1);
[coefV,errorEst]=polyfit(xRej(:),yRej(:),1);
tmp = polyparci(coefV,errorEst,0.95); % gives lower and upper CI
ci95 = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)]; % only keep the diff
tmp = polyparci(coefV,errorEst,0.99);
ci99 = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)];



%%%%% is there a difference between static and motion?
% compare the slope for motion and static (only 3 datapoints)
clear data y;
y(:,:,1) = f1Amp([11 13 15],:)'; % 2.5 Hz static
y(:,:,2) = f1Amp([16 21 20],:)'; % 2.5 Hz moving
y(:,:,3) = f1Amp([7 8 9],:)'; % 5 Hz static
y(:,:,4) = f1Amp([17 22 19],:)'; % 5 Hz moving

for cond=1:4
    data(cond,:) = reshape(y(:,:,cond),1,length(keepSbj)*3);
end
for cond=1:4
    if cond<3
        x = repmat([12.5 50 87.5],length(keepSbj),1);
    else
        x = repmat([25 50 75],length(keepSbj),1);
    end
    [coefV,errorEst]=polyfit(x(:)',data(cond,:),1);
    coef(:,cond)=coefV;
    tmp = polyparci(coefV,errorEst,0.95); % gives lower and upper CI
    ci95(:,cond) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)]; % only keep the diff
    tmp = polyparci(coefV,errorEst,0.99); 
    ci99(:,cond) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)];
end






%%%%%%%%%%% behaviour
clear x y data coef ci95 ci99;
load('fullRatings9.mat')

%%%%% compare average static vs moving ratings
load fullRatings9.mat
flicker = squeeze(mean(mean(tabStatic,1),2));
for ss=1:9
    move(ss) = mean(nonzeros(tabMot(:,:,ss)));
end
moving = moving';
[P,H,STATS] = signrank(flicker,move); % signrank Wilcoxon signed rank test 


static = mean(tabStatic,3); moving = mean(tabMot,3);

%%%%% Effect of DC on ratings for flicker?
clear coef;
x = repmat(dcVal,length(tabStatic),1);
y = static(1,1:5,:)'; % 10hz static
data(1,:) = y(:);
y = static(2,1:5,:)'; % 5Hz static
data(2,:) = y(:);
y = static(3,1:5,:)'; % 2.5Hz static
data(3,:) = y(:);

fqCond = {'10','5','2.5'};
for freq = 1: 3
    y = squeeze(tabStatic(freq,1:5,:))'; 
    [coefV,errorEst]=polyfit(x(:),y(:),1);
    coef(:,freq)=coefV;
    tmp = polyparci(coefV,errorEst,0.95); % gives lower and upper CI
    ci95(:,freq) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)]; % only keep the diff
    tmp = polyparci(coefV,errorEst,0.99); 
    ci99(:,freq) = [tmp(2,1) - coefV(1) tmp(2,2) - coefV(1,2)];
    result = [fqCond{freq}  'slope ' num2str(coef(1,freq)) ' +/-' num2str(ci95(1,freq)) ' +/-' num2str(ci99(1,freq))];
    disp(result)
    result = [fqCond{freq}  'intercept ' num2str(coef(2,freq)) ' +/-' num2str(ci95(2,freq)) ' +/-' num2str(ci99(2,freq))];
    disp(result)
end

%%%%% Effect of DC on ratings for moving stim?
clear coef yM xM tmp;
xM = repmat(dcVal(2:4),length(tabMot),1);
yM = squeeze(tabMot(2,2:4,:))'; % 5Hz
[coef,error]=polyfit(xM,yM,1);
tmp = polyparci(coef,error,0.95);
[ coef ; tmp(2,1) - coef(1) tmp(2,2) - coef(1,2)] % 1st line is coef and intercept, 2nd is 95CI
clear coef yM xM tmp;
xM = repmat(dcVal([1 3 5]),length(tabMot),1);
yM = squeeze(tabMot(3,[1 3 5],:))'; % 5Hz
[coef,error]=polyfit(xM,yM,1);
tmp = polyparci(coef,error,0.95);
[ coef ; tmp(2,1) - coef(1) tmp(2,2) - coef(1,2)] % 1st line is coef and intercept, 2nd is 95CI
