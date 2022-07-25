clearvars

% load models
load('eccModel.mat');

% remove 0 in the model matrix
[nbElec, tot] = size(eccModelRef.amp);
tmpAmp = nonzeros(getfield(eccModelRef,'amp'));
tmpAmpReshape = reshape(tmpAmp,nbElec-1,length(tmpAmp)/(nbElec-1));
tmpCos = nonzeros(getfield(eccModelRef,'cos'));
tmpCosReshape = reshape(tmpCos,nbElec-1,length(tmpCos)/(nbElec-1));

minModel.amp = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.amp(2:end,:) = tmpAmpReshape;
minModel.cos = zeros(nbElec,length(tmpAmp)/(nbElec-1));
minModel.cos(2:end,:) = tmpCosReshape;
% amplitude = abs(cos), no need to get it from matrix
% also there is no phase so sin = zeros
minModel.sin = zeros(size(minModel.cos));

%Normalize column power
%Note! this uses implicit binary singleton expansion
minModel.cos = minModel.cos./sum(minModel.cos.^2);

dataDir = '/Volumes/Amrutam/Marlene/JUSTIN/DutyCycle/data/Axx/';

listData = dir([dataDir '*.mat']);
cfg.layout = 'biosemi128.lay';
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
dcVal = [12.5 25 50 75 87.5];

keepSbj = [1:5 7:11 13:15 17:18 20];
% reject S7 and S13 S17, S20 same as S03
% S1=pilote so reject numS-1 (6 12 16 19)

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
   
    noiseData = [];
    sigData = [];
    for cond = 1:length(Axx)
        sigData = [sigData Axx(cond).sin(:,Axx(cond).harmIdx) Axx(cond).cos(:,Axx(cond).harmIdx)];
        noiseData = [noiseData Axx(cond).sin(:,Axx(cond).harmIdx-1) Axx(cond).cos(:,Axx(cond).harmIdx-1)...
            Axx(cond).sin(:,Axx(cond).harmIdx+1) Axx(cond).cos(:,Axx(cond).harmIdx+1)];
    end
    noiseCov = cov(noiseData');
    sigCov = cov(sigData');
    figure; set(gcf, 'Position', [0 0 1000 500])
    subplot(1,2,1);imagesc(sigCov);colorbar;title(['signal covariance S' num2str(ss)]);
    subplot(1,2,2);imagesc(noiseCov);colorbar;title('noise covariance');
    saveas(gcf,['figures' filesep 'covS' num2str(ss)],'png')
end

