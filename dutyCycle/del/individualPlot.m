
% how do I do the regression? on average or on each sbj separately?
% how do I compare to noise? using average? each sbj to its own noise level?

clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork  
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions
ft_defaults
dataDir = '/Users/marleneponcet/Documents/data/dutyCycle/Axx/';


listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};

keepSbj = [1:5 7:11 13:15];

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    dataSbj(:,ss) = Axx;
end


pickElec = 23;

%%% get just the 1f1
for elec=1:length(pickElec)

    for ss=1:length(keepSbj)
        for cond=1:15
            fAmp(cond,ss)=dataSbj(cond,ss).amp(pickElec(elec),(dataSbj(cond,ss).i1f1));
        end
        for cond=16:22 % for motion take the 2f1
            fAmp(cond,ss)=dataSbj(cond,ss).amp(pickElec(elec),(dataSbj(cond,ss).i1f1*2-1));
        end
    end

end


%%% plot each participant
figure;hold on;
for ss=1:13
    subplot(3,5,ss); hold on
    for freq=1:3
        plot(fAmp((freq-1)*5+1:freq*5,ss),'Linewidth',2)
    end
    plot(3,fAmp(18,ss),'^:','Linewidth',1)
    plot([2 3 4],fAmp([17 22 19],ss),'^:','Linewidth',1)
    plot([1 3 5],fAmp([16 21 20],ss),'^:','Linewidth',1)
end

dcVal = [12.5 25 50 75 87.5];

%%% plot all participants for each cond separately
figure;hold on;
for freq=1:3
    subplot(2,2,freq); hold on
    plot(dcVal,fAmp((freq-1)*5+1:freq*5,:),'Linewidth',2)
end







