%%% normalise the Fourier spectrum depending on a simulated square function
%%% if you want to normalise for an impulse function, would need to set the
%%% amplitude of the impulse for a max = 1 (numbers are very low otherwise)



clearvars; 
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated

load('16sbjDC.mat')
load('simulDC.mat')
listFq = avData(1).freq;
% maxFreq = find(listFq==85/8*4);
f1list = [avData(1).i1f1 avData(6).i1f1 avData(11).i1f1];
dcVal = [12.5 25 50 75 87.5];
col={'b','r','g'};

% determine all the harmonics
for fq = 1:3
    harm(fq).fqIndex = determineFilterIndices( 'nf1low49', avData(1+5*(fq-1)).freq, avData(1+5*(fq-1)).i1f1);
end

% simulDC are stil imaginary nb, change to amplitudes first
squareAmp = abs(squareFFT); 

% just work with Oz freq
dataOz = zeros(size(squareFFT));
for fq = 1:3
    for dc = 1:5
        dataOz(:,fq,dc) = avData(dc+(fq-1)*5).amp(23,:);
    end
end

% sum all the freq harm < maxFreq
for fq = 1:3
    for dc = 1:5
        sumData(fq,dc) = sum(dataOz(harm(fq).fqIndex,fq,dc));
        sumSquare(fq,dc) = sum(squareAmp(harm(fq).fqIndex,fq,dc));
    end
end
normSqData = sumData ./ sumSquare;



figure; hold on;
% plot f1 amplitude across conditions 
subplot(2,2,1); hold on;
for fq = 1:3
    plot(dcVal,squeeze(dataOz(f1list(fq),fq,:)),['.-' col{fq}],'MarkerSize',40,'Linewidth',2);
end
xlim([0 100]);
% xticks([0:12.5:100]);
legend('10','5','2.5','Location','Best')
xlabel('Duty Cycle')
ylabel('SSVEP amplitude')
title('1f1')
% plot sum all harm across conditions 
subplot(2,2,2); hold on;
for fq = 1:3
    plot(dcVal,sumData(fq,:),['.-' col{fq}],'MarkerSize',40,'Linewidth',2);
end
xlim([0 100]);
title('sum all harm < 50Hz')
% plot all harm amp / square
subplot(2,2,3); hold on;
for fq = 1:3
    plot(dcVal,sumSquare(fq,:),['.-' col{fq}],'MarkerSize',40,'Linewidth',2);
end
xlim([0 100]);
title('sum square')
% plot all harm amp / square
subplot(2,2,4); hold on;
for fq = 1:3
    plot(dcVal,normSqData(fq,:),['.-' col{fq}],'MarkerSize',40,'Linewidth',2);
end
xlim([0 100]);
title('all harm / square')




