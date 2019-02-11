%%%%%% test stuff

%%%%%%%%%%%%
%% RMS
for ss=1:3
    for chan=1:2
        for tt=1:10
            test(ss,chan,tt) = randi(10);
        end
    end
end
% same as:
test = randi([0, 10], [3,2,10]); % subj, chan, time
% sbj average then rms time
avTest = squeeze(mean(test));
rms(avTest,2)

% rms time then average sbj
timeTest = rms(test,3);
mean(timeTest)

%%%%%%%%%%%%%%%
%% test filtering neighboors
clearvars
load('diff.mat')
filtIdx = determineFilterIndices( 'nf1low50', diff(1,1).freq, diff(1,1).i1f1 );
filtVoisin = [filtIdx-1 filtIdx+1 ];

%Create a logical matrix selecting frequency components.
filtMat = false(size(diff(1,1).amp));
filtMat(:,filtVoisin) = true;

%Combine the filter and sig vaules with a logical AND.
% condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
diff(1,1).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
cfg.activeFreq =  diff(1,1).activeFreq;

[ diff(1,1).baselineWave ] = filterSteadyState( cfg, diff(1,1) );


%%%%%%%%%%%
%% regress
% multiple linear regression
lin = [1 2 1]';
temp = [1 0 1]';
spat = [-3 2 1]';
actSig = 0.5*lin + 3*temp + 2*spat ;
regCoef = [ones(size(lin)) lin temp spat spat.*temp]; % as in a+bLin+cTemp+dSpat+eTemp*Spat
regress(actSig,regCoef) % Y,X !!

load('results.mat') % average
for elec=1:size(pred(1,1).filteredWave,1)
    longAM = pred(1,1).filteredWave(elec,:)';
    regCoef = [pred(1,2).filteredWave(elec,:)' pred(1,1).filteredWave(elec,:)'-pred(1,3).filteredWave(elec,:)' pred(1,1).filteredWave(elec,:)'-pred(1,4).filteredWave(elec,:)']; % as linear + spatial + temporal
    longCoef(elec,:) = regress(longAM,regCoef); 
    shortAM = pred(1,6).filteredWave(elec,:)';
    regCoef = [pred(1,7).filteredWave(elec,:)' pred(1,1).filteredWave(elec,:)'-pred(1,8).filteredWave(elec,:)' pred(1,1).filteredWave(elec,:)'-pred(1,9).filteredWave(elec,:)']; % as linear + spatial + temporal
    shortCoef(elec,:) = regress(shortAM,regCoef); 
end
figure;
for fact=1:3
    subplot(3,1,fact)
    plotTopo(longCoef(:,fact),cfg.layout)
    colorbar
end
figure;
for fact=1:3
    subplot(3,1,fact)
    plotTopo(shortCoef(:,fact),cfg.layout)
    colorbar
end
for numCond=1:4 % long/short E1 E2
for elec=2:size(pred(1,1).filteredWave,1) % first electrode is the reference 0
    longAM = pred(1,1+5*(numCond-1)).filteredWave(elec,:)';
    regCoef = [pred(1,2+5*(numCond-1)).filteredWave(elec,:)' pred(1,1+5*(numCond-1)).filteredWave(elec,:)'-pred(1,3+5*(numCond-1)).filteredWave(elec,:)' pred(1,1+5*(numCond-1)).filteredWave(elec,:)'-pred(1,4+5*(numCond-1)).filteredWave(elec,:)']; % as linear + spatial + temporal
    coefFact(numCond,elec,:) = regress(longAM,regCoef); 
end
end
% reconstruct signal from the coef

plot(longCoef(elec,2)*(pred(1,16).filteredWave(elec,:)'-pred(1,18).filteredWave(elec,:)'))
plot(longCoef(elec,1)*pred(1,17).filteredWave(elec,:)' + longCoef(elec,2)*(pred(1,16).filteredWave(elec,:)'-pred(1,18).filteredWave(elec,:)') + longCoef(elec,3)*(pred(1,16).filteredWave(elec,:)'-pred(1,19).filteredWave(elec,:)'))
hold on 
plot(pred(1,16).filteredWave(elec,:)')
plot(pred(1,18).filteredWave(elec,:)')

% multivariate linear regression
lin = [1 2 1; 3 5 1]';
temp = [1 0 1; 3 0 1]';
spat = [-3 2 1; 0 3 3]';
actSig = 0.5*lin + 3*temp + 2*spat ;
regCoef = {ones(size(lin)) lin temp spat spat.*temp}; % as in a+bLin+cTemp+dSpat+eTemp*Spat
mvregress(regCoef,actSig) % X, Y!! 
[beta,Sigma,E,CovB,logL] = mvregress(X,Y);

clear regCoef
regCoef(1,:,:) = temp;
regCoef(2,:,:) = spat;
mvregress(regCoef,actSig)



