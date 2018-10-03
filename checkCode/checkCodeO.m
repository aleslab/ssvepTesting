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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% test Axx
load('AxxTestS4')
AxxToShift = Axx(2);
AxxNoChange = Axx(3);

% compute the new Axx with a shift in time of 200 ms (half cycle)
% first for the wave
shift.wave = circshift(AxxToShift.wave,[0 length(AxxToShift.wave)/2]);
% then for the cos and sin
% for this needs the phase shift calculated for each frequency
% phase shift = shift (200 ms) / periods * 360 degrees
periods = 1 ./ AxxToShift.freq;
% phaseShift = 0.2 ./ periods .* 360;
phaseShift = AxxToShift.time(end/2+1)/1000 ./ periods .* 360; 

% trigonometry.. 
% cos' = A*cos(theta + phi) = cos(theta)*alpha*cos(phi)
% (theta = angle difference, phi = original angle)
% cos' = cosd(phaseShift) * orignialCos - sind(phaseShift) * orignialSin
% Similar stuff for sin' but different "identification"
% sin' = sind(phaseShift) * originalCos + cosd(phaseShift) * originalSin
shift.sin = sind(phaseShift) .* AxxToShift.cos + cosd(phaseShift) .* AxxToShift.sin;
shift.cos = cosd(phaseShift) .* AxxToShift.cos - sind(phaseShift) .* AxxToShift.sin;

% now do the addition
motion.wave = AxxNoChange.wave +shift.wave;
motion.sin = AxxNoChange.sin + shift.sin;
motion.cos = AxxNoChange.cos + shift.cos;
motion.amp = sqrt(motion.sin.^2 + motion.cos.^2);

iChan=23;

% construct wave average from sin/cos for chan 23
dft = dftmtx(AxxToShift.ndft);
dft = dft(:,1:AxxToShift.nfr);
waveRecon = real(dft(:,:))*motion.cos(iChan,:)' - imag(dft(:,:))*motion.sin(iChan,:)';
waveRecon=reshape(waveRecon,AxxToShift.nt,[])';
motion.waveRecon(iChan,:)=mean(waveRecon,1);

figure; hold on;
plot(motion.wave(iChan,:),'r')
plot(AxxNoChange.wave(iChan,:),'b')
plot(shift.wave(iChan,:),'g')
plot(motion.waveRecon(iChan,:),'k')
legend('A+Bshift','A','Bshift','recon')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('testAv') % this is 4 sbj of LRshortDC sbjprediction with 2 conditions
load('testAllCond') % this is 4 sbj of LRshortDC sbjprediction with all 12 conditions
testAv = testAllCond;

%%%%%%%%%%%% test ind sbj
for ss=1:4
    figure; hold on;
for cond=1:12
subplot(3,4,cond);hold on
plot(testAv(ss,cond).data.wave(23,:));
iChan=23;
dft = dftmtx(testAv(ss,cond).data.ndft);
dft = dft(:,1:testAv(ss,cond).data.nfr);
waveRecon = real(dft(:,:))*testAv(ss,cond).data.cos(iChan,:)' - imag(dft(:,:))*testAv(ss,cond).data.sin(iChan,:)';
waveRecon=reshape(waveRecon,testAv(ss,cond).data.nt,[])';
wave=mean(waveRecon,1);
plot(wave)
legend({'wave1cond','reconstruct'})
end
end

%%%%%%%%%%%% test average
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
    clear waveRecon;
    dft = dftmtx(testAv(1,cond).data.ndft);
    dft = dft(:,1:testAv(1,cond).data.nfr);
    waveRecon = real(dft(:,:))*groupAv(cond).cos(iChan,:)' - imag(dft(:,:))*groupAv(cond).sin(iChan,:)';
    waveRecon=reshape(waveRecon,testAv(cond).data.nt,[])';
    groupAv(cond).waveRecon(iChan,:)=mean(waveRecon,1);
end

% plot
for cond=1:2
    figure;hold on;
plot(groupAv(cond).wave(iChan,:),'b','LineWidth',2)
plot(groupAv(cond).waveRecon(iChan,:),'r','LineWidth',2)
for ff=1:size(testAv,1)
    plot(dataWave(iChan,:,cond,ff))
end
legend('waveAverage','reconstructWave','each sbj')
end












%%%%%%%%%%%% test filter
dataIn = {'/Users/marleneponcet/Documents/data/LongRangeV2/', '/Users/marleneponcet/Documents/data/LRshortDC/V2/'};
labelIn = {'75DC','25DC'};
colDiff = {'r','g','b','c'};
elec = [23 19]; % 23=Oz, 19=Pz, 3=CPz, A7/B4-36=P3/P4, FCz=87, C3=115, C4=54
saveplot = 0;

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