% Compute average of Axx with Bootstrap + power
% sum of power computed on sum of amplitude harmonics (for the group)
% Different from sum of power per sbj then average
% Harmonics NOT normalised yet

clearvars;
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults


% load individual data
dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';
listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
for ss = 1:length(listData)
    clear Axx;
    load([dataDir listData(ss).name]);
    dataSbj(:,ss) = Axx;
end
[avData, proj_Amp] = averageAxxWithStd(dataSbj);

% compute sum of power for signal and noise
for cond=1:length(dataSbj)
    avData(cond).sumPower = sum(avData(cond).amp(:,avData(cond).harmIdx).^2,2);
    avData(cond).sumPowerNoise = sum(avData(cond).amp(:,[avData(cond).harmIdx-1 avData(cond).harmIdx+1]).^2,2) /2;
end
% bootstrap for getting confidence intervals
nbBoot = 1000;
for cond=1:length(dataSbj)
    power = zeros(nbBoot,avData(1).nchan);powerNoise = zeros(nbBoot,avData(1).nchan);
    for bb=1:nbBoot
        if mod(bb,100)==0
            fprintf('condition%d boostrap%d \n' ,cond,bb)
        end
        % pick randomly 16 sbj with replacement
        pickSS = randi(length(listData),1,length(listData));
        tmpSin = zeros(length(listData),128,577); tmpCos = zeros(length(listData),128,577);
        for ss=1:length(listData)
            tmpSin(ss,:,:) = dataSbj(cond,pickSS(ss)).sin;
            tmpCos(ss,:,:) = dataSbj(cond,pickSS(ss)).cos;
        end
        % compute power for all freq
        powerAll = squeeze(mean(tmpSin)).^2 + squeeze(mean(tmpCos)).^2;
        % sum only the harmonics
        power(bb,:) = sum(powerAll(:,avData(cond).harmIdx),2);
        powerNoise(bb,:) = sum(powerAll(:,[avData(cond).harmIdx-1 avData(cond).harmIdx+1]),2) /2;
    end
    avData(cond).ciPower = prctile(power,[2.5 97.5])'; 
    avData(cond).ciPowerNoise = prctile(powerNoise,[2.5 97.5])';
    avData(cond).stdPower = std(power)'; 
    avData(cond).stdPowerNoise = std(powerNoise)';
    avData(cond).stdAmp = std(sqrt(power))';
    avData(cond).stdAmpNoise = std(sqrt(powerNoise))';
    avData(cond).ciAmp = prctile(sqrt(power),[2.5 97.5])'; 
    avData(cond).ciAmpNoise = prctile(sqrt(powerNoise),[2.5 97.5])';
    % note that std(power)/sumSq = std(power/sumSq) so normalisation can be
    % done later on the stdPower
end

% % normalisation term
% % load the amplitudes for a square function
% load('squareFFT.mat')
% for cond=1:length(avData)
%     sumSq(cond) = sum(sqAllFFT(avData(cond).harmIdx,cond).^2);
% end

save('dataDCavRef','avData','proj_Amp','cfg')
