clearvars;

%% create stimuli
% moving stim similar to the exp but 4 locations so that it is going to the
% right (net motion different from 0)

%%% 2.5 Hz
f1 = [ones(40,50) zeros(40,150)];
mouv(:,:,3) = [f1; circshift(f1,50,2); circshift(f1,100,2); circshift(f1,150,2)];
figure;imagesc(mouv(:,:,3))

%%% 5 Hz
f1 = [ones(20,50) zeros(20,150)];
cycle4 = [f1; circshift(f1,50,2); circshift(f1,100,2); circshift(f1,150,2)];
mouv(:,:,2) = repmat(cycle4,2,1);
figure;imagesc(mouv(:,:,2))


%%% 10 Hz
f1 = [ones(10,50) zeros(10,150)];
cycle4 = [f1; circshift(f1,50,2); circshift(f1,100,2); circshift(f1,150,2)];
mouv(:,:,1) = repmat(cycle4,4,1);
figure;imagesc(mouv(:,:,1))

fq=[10 5 2];

for fff = 1:size(mouv,3)
    stim = mouv(:,:,fff);
[total_energy, motion_energy, energy_flicker, motion_contrast, peak] = energy(stim);

            flickerMouv(fff) = energy_flicker;
            motionMouv(fff) = motion_energy;
            
% Plot the stimulus
figure; 
subplot(2,1,1);
imagesc(stim); 
colormap(gray);
axis off
caxis([0 1.0]);
axis equal
title('Stimulus');

subplot(2,1,2);
imagesc(motion_contrast); 
colormap(gray);
axis off
caxis([-peak peak]); colorbar;
axis equal
title('Normalised Motion Energy');

saveas(gcf,['realMotion' num2str(fq(fff)) 'Hz'],'png' );
end

save('sumMouv.mat','flickerMouv','motionMouv')

