close all;
clearvars;


addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFuntions
ft_defaults

% addpath C:\Users\Marlene\Documents\git\fieldtrip-aleslab-fork
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\svndlCopy
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\biosemiUpdated
% addpath C:\Users\Marlene\Documents\git\ssvepTesting\commonFuntions
% ft_defaults
% dataIn = {'C:\Users\Marlene\Documents\git\dataLR\LRlongDC\', 'C:\Users\Marlene\Documents\git\dataLR\LRshortDC\V2\'};


dataIn = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};
labelIn = {'75DC','25DC'};
colDiff = {'r','g','b','c'};
elec = [23 19]; % 23=Oz, 19=Pz, 3=CPz, A7/B4-36=P3/P4, FCz=87, C3=115, C4=54
saveplot = 0;


%%%%%%%%%%%% test 1 sbj
figure; hold on;
plot(testAv(1,1).data.wave(23,:));
iChan=23;
dft = dftmtx(testAv(1,1).data.ndft);
dft = dft(:,1:testAv(1,1).data.nfr);
waveRecon = real(dft(:,:))*testAv(1,1).data.cos(iChan,:)' - imag(dft(:,:))*testAv(1,1).data.sin(iChan,:)';
waveRecon=reshape(waveRecon,testAv(1,1).data.nt,[])';
wave=mean(waveRecon,1);
plot(wave)
legend({'wave1cond','reconstruct'})


%%%%%%%%%%%% test average

load('testAv')

% reshape data
for ff=1:size(testAv,1)
    for cond=1:size(testAv,2)
        dataWave(:,:,cond,ff) = testAv(ff,cond).data.wave;
        dataAmp(:,:,cond,ff) = testAv(ff,cond).data.amp;
        dataSin(:,:,cond,ff) = testAv(ff,cond).data.sin;
        dataCos(:,:,cond,ff) = testAv(ff,cond).data.cos;
    end
end

% average
for cond=1:2
groupAv(cond).wave = mean(dataWave(:,:,cond,:),4);
groupAv(cond).sin = mean(dataSin(:,:,cond,:),4);
groupAv(cond).cos = mean(dataCos(:,:,cond,:),4);
groupAv(cond).amp = sqrt(groupAv(cond).sin.^2 + groupAv(cond).cos.^2);
end

% construct wave average from sin/cos for chan 23
iChan=23;
for cond=1:2
    dft = dftmtx(testAv(cond).data.ndft);
    dft = dft(:,1:testAv(cond).data.nfr);
    waveRecon = real(dft(:,:))*testAv(cond).data.cos(iChan,:)' - imag(dft(:,:))*testAv(cond).data.sin(iChan,:)';
    waveRecon=reshape(waveRecon,testAv(cond).data.nt,[])';
    wave=mean(waveRecon,1);
end

% plot
cond=1;figure;hold on;
plot(groupAv(cond).wave(iChan,:),'b','LineWidth',2)
plot(wave,'r','LineWidth',2)
for ff=1:size(testAv,1)
    plot(dataWave(iChan,:,cond,ff))
end

legend('waveAverage','reconstructWave','each sbj')


%%%%%%%%%%%% test filter

for dd=1:length(dataIn) % 2 experiments
    clear avPredictions sbj diff avInteractions;
    %%% load the sbj predictions
    load([dataIn{dd} 'sbjprediction.mat'])
    load([dataIn{dd} 'NLinteraction.mat'])
    
    % do the average
    avPredictions = averageSbj(sbj);
    avInteractions = averageAxx(interaction');
    
    %%%%%% compute differences
    pickCond = [2:5 ; 7:10];
    for ss=1:size(sbj,1)
        ppp=0;
        for cc=1:size(pickCond,1)
            condCompare = pickCond(cc,:);
            fixCond = condCompare(1) - 1;
            for cond=1:length(condCompare)
                ppp=ppp+1;
                test(ss,ppp).data = computeDiff(sbj(ss,fixCond).data,sbj(ss,condCompare(cond)).data);
                diff(cc,cond) = computeDiff(avPredictions(fixCond),avPredictions(condCompare(cond)));
            end
        end
    end
    figure;hold on;
    for ss=1:7
        plot(test(ss,1).data.wave(23,:))
    end
    avTest = averageSbj(test);
    
    plot(avTest(1).wave(23,:),'LineWidth',3)
    plot(diff(1,1).wave(23,:),'LineWidth',3)
    
    figure;hold on;
    for cc=5:8
        plot(avTest(cc).wave(23,:))
        plot(diff(2,cc-4).wave(23,:))
    end
    ylim([-3 3])
    
    
    % warning about filter is NORMAL
    for condIdx=1:length(avTest)
        filtIdx = determineFilterIndices( 'low49', avTest(condIdx).freq, avTest(condIdx).i1f1 );
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(avTest(condIdx).amp));
        filtMat(:,filtIdx) = true;
        
        %Combine the filter and sig vaules with a logical AND.
        % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        avTest(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
        cfg.activeFreq =  avTest(condIdx).activeFreq;
        
        [ avTest(condIdx).filteredWave ] = filterSteadyState( cfg, avTest(condIdx) );
    end
    
    figure;hold on;
    pp=0;
    for cc=[5 8]
        pp = pp+1;
        col={'r','b'};
        plot(avTest(cc).time,avTest(cc).filteredWave(23,:),col{pp},'LineWidth',2)
        plot(avTest(cc).time,avTest(cc).wave(23,:),col{pp})
    end
    legend({'filtered','unfiltered'})
    
    
    
    %%%%%%%%% CHECK interaction
    figure;hold on;
    for ss=1:7
        plot(interaction(ss,1).wave(23,:))
    end
    plot(avInteractions(1).wave(23,:),'LineWidth',3)
    figure;hold on;
    for cc=4:6
        plot(avInteractions(cc).wave(23,:))
    end
    
    for condIdx=1:length(avInteractions)
        filtIdx = determineFilterIndices( 'low49', avInteractions(condIdx).freq, avInteractions(condIdx).i1f1 );
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(avInteractions(condIdx).amp));
        filtMat(:,filtIdx) = true;
        
        %Combine the filter and sig vaules with a logical AND.
        % condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        avInteractions(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
        cfg.activeFreq =  avInteractions(condIdx).activeFreq;
        
        [ avInteractions(condIdx).filteredWave ] = filterSteadyState( cfg, avInteractions(condIdx) );
    end
    
    figure;hold on;
    pp=0;
    for cc=[1 3]
        pp = pp+1;
        col={'r','b'};
        plot(avInteractions(cc).time,avInteractions(cc).filteredWave(23,:),col{pp},'LineWidth',2)
        plot(avInteractions(cc).time,avInteractions(cc).wave(23,:),col{pp})
    end
    legend({'filtered','unfiltered'})
end