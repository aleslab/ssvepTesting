%% lod data 
load('testAv')
 
%% analyze
% reshape data
for ff=1:size(testAv,1)
    for cond=1:size(testAv,2)
        dataWave(:,:,cond,ff) = testAv(ff,cond).data.wave;
        dataAmp(:,:,cond,ff) = testAv(ff,cond).data.amp;
        dataSin(:,:,cond,ff) = testAv(ff,cond).data.sin;
        dataCos(:,:,cond,ff) = testAv(ff,cond).data.cos;
        
        %Recon the wave from sin/cos in each subj
        dft = dftmtx(testAv(ff,cond).data.ndft);
        dft = dft(:,1:testAv(ff,cond).data.nfr);
        
        for iChan = 1:testAv(ff,cond).data.nchan,                    
            tmpWave = real(dft(:,:))*testAv(ff,cond).data.cos(iChan,:)' - imag(dft(:,:))*testAv(ff,cond).data.sin(iChan,:)';
            tmpWave=reshape(tmpWave,testAv(ff,cond).data.nt,[])';
            dataWaveRecon(iChan,:,cond,ff) = mean(tmpWave,1);
        end
        
    end
end
 
% average
for cond=1:2
groupAv(cond).wave = mean(dataWave(:,:,cond,:),4);
groupAv(cond).sin = mean(dataSin(:,:,cond,:),4);
groupAv(cond).cos = mean(dataCos(:,:,cond,:),4);
groupAv(cond).amp = sqrt(groupAv(cond).sin.^2 + groupAv(cond).cos.^2);
end
 
% construct wave average from sin/cos for all channels
for cond=1:2
    dft = dftmtx(testAv(1,cond).data.ndft);
    dft = dft(:,1:testAv(1,cond).data.nfr);
    
    for iChan = 1:testAv(1,cond).data.nchan
 
        tmpRecon = real(dft(:,:))*groupAv(cond).cos(iChan,:)' - imag(dft(:,:))*groupAv(cond).sin(iChan,:)';
        tmpRecon=reshape(tmpRecon,testAv(1,cond).data.nt,[])';
        groupAv(cond).waveRecon(iChan,:)=mean(tmpRecon,1);
    end
    
end
 
%% Plot
% plot
iChan=23;
    
for cond = [1 2];
figure(100+cond);clf;
plot(groupAv(cond).waveRecon(iChan,:),'r','LineWidth',2)
title(['Group Average participants for condition: ' num2str(cond) ])
hold on;
plot(groupAv(cond).wave(iChan,:),'b','LineWidth',2)
for ff=1:size(testAv,1)
    plot(dataWave(iChan,:,cond,ff))
end
 
legend('waveRecon','groupAverage','each sbj')
 
figure(200+cond);
clf;
for iP = 1:size(dataWave,4),
subplot(size(dataWave,4),1,iP);
title(['Individual participants for condition: ' num2str(cond) ])
plot(dataWave(iChan,:,cond,iP),'k');
hold on
plot(dataWaveRecon(iChan,:,cond,iP),'--k')
plot(dataWave(iChan,:,cond,iP)-dataWaveRecon(iChan,:,cond,iP),'-b')
end
legend('data in Wave', 'Data from sin/cos','Difference');
 
 
 
end
 
 