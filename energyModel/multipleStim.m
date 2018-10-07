
% add other durations: 100, 200, 400 for around 2.5,5,10Hz
% this should be changed in the timeON loop (+ repmat, etc.)

%% create stimuli
% should be 200 x 160 
% 200 = 8 degrees, separation of less than a degree = 8+8 = 16 = 0.64
% Will have 2 stim each of 46 = 1.84 degrees
% 160 = 1.6 s, each step = 10 ms, half cycle (one side) = 20, full cycle (2
% sides = motion) = 40

clearvars;


stimWidth = 46;
totWidth = 200;
midSep = 8;

cycleLength = [20 40 80]; % 200 400 800ms
duration = 160;
for cc=1:length(cycleLength)
    timeON(cc,:) = [10:10:100]/100 * cycleLength(cc)/2;
end

% location static stimulus
locS1(1,:) = totWidth/2-stimWidth/2+1: totWidth/2+stimWidth/2; locS2(1,:) = locS1(1,:);
% location short-range
locS1(2,:) = stimWidth+1:totWidth/2-midSep; locS2(2,:) = totWidth/2+midSep+1:totWidth-stimWidth;
% location long-range
locS1(3,:) = 1:stimWidth; locS2(3,:) = totWidth-stimWidth+1:totWidth;
for tt=1:size(locS1,1)
    for cc = 1:length(cycleLength)
        for on=1:length(timeON)
            oneCycle = zeros(cycleLength(cc),totWidth);
            oneCycle(1:timeON(cc,on),locS1(tt,:)) = ones(timeON(cc,on),stimWidth);
            oneCycle(cycleLength(cc)/2+1:cycleLength(cc)/2+timeON(cc,on),locS2(tt,:)) = ones(timeON(cc,on),stimWidth);
            allStim(:,:,cc,tt,on) = repmat(oneCycle,duration/cycleLength(cc),1);
        end
        figure;imagesc(allStim(:,:,cc,tt,on))
    end
end





%% add loop here for all the stim to test
for tt = 1:size(allStim,4)
    for cc=1:size(allStim,3)
        for ss=1:size(allStim,5)
            stim = allStim(:,:,cc,tt,ss);
            
            [total_energy, motion_energy, energy_flicker, motion_contrast, peak] = energy(stim);
            
            flickerEnergy(ss,cc,tt) = energy_flicker;
            motionEnergy(ss,cc,tt) = motion_energy;
            
            % if ss == 1 || ss==size(allStim,3)/2 || ss==size(allStim,3)
            if ss==size(allStim,5)/2
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
                
                % saveas(gcf,['stim' num2str(tt) 'dc' num2str(ss/size(allStim,3)*100) 'png'],'png' )
                saveas(gcf,['stim' num2str(tt) 'Fq' num2str(cc) ],'png' )
            end
            
            
        end
    end
end

load('sumMouv.mat')
fq = [10 5 2.5];

figure; 
for cc=1:size(allStim,3)
    subplot(1,3,cc);hold on;
    for tt = 1:size(allStim,4)
        plot(flickerEnergy(:,cc,tt))
    end
    plot(10,flickerMouv(cc),'Marker','*')
    title([num2str(fq(cc)) 'Hz'])
    if cc ==1; ylabel('flickerEnergy (= average total energy)'); end
end
legend('static','SR','LR','motion')
saveas(gcf,'Flicker summary','png' )

figure;
for cc=1:size(allStim,3)
    subplot(1,3,cc);hold on;
    for tt = 1:size(allStim,4)
        plot(motionEnergy(:,cc,tt))
    end
    plot(10,motionMouv(cc),'Marker','*')
    title([num2str(fq(cc)) 'Hz'])
    if cc ==1; ylabel('motionEnergy'); end
end
legend('static','SR','LR','motion')

saveas(gcf,'Motion summary','png' )

