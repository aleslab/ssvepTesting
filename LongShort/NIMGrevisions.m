
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check nb of rejected trials per condition
clearvars;
dataPath = {'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E1LRshortDC\Axx\',...
    'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E2LRlongDC\Axx\'};
% max nb of epochs = 72 in shortDC (18blk*4epochs)
% = 102 in longDC (17blk * 6 epochs)
maxEpoch = [72, 102];

for exp=1:length(dataPath)
    clear totTrial rejTrial
    axxFiles = dir([dataPath{exp} '*.mat']);
    for sbj=1:length(axxFiles)
        clear Axx;
        load([dataPath{exp} axxFiles(sbj).name]);
        for cond=1:length(Axx)
            totTrial(sbj,cond) = length(Axx(cond).cfg.trials); % this is just a check
            rejTrial(sbj,cond) =( maxEpoch(exp) - length(Axx(cond).cfg.trials)) / maxEpoch(exp)*100;
        end
    end
    clear g rejTrial2;
%     figure;
%     rejTrial2 = rejTrial(:,1:6)';
%     g(1,1) = gramm('x',repmat(1:6,1,length(rejTrial)),'y',rejTrial2(:),'ymin',repmat(0,length(rejTrial2(:)),1),'ymax',repmat(60,length(rejTrial2(:)),1));
%     g(1,1).stat_boxplot();
%     rejTrial3 = rejTrial(:,7:12)';
%     g(2,1) = gramm('x',repmat(1:6,1,length(rejTrial)),'y',rejTrial3(:),'ymin',repmat(0,length(rejTrial3(:)),1),'ymax',repmat(60,length(rejTrial3(:)),1));
%     g(2,1).stat_boxplot();
%     g.draw();  
    figure;
    rejTrial2 = rejTrial(:,1:6)'; % {'AM','sL','sR','sim','dL','dR'}
    g = gramm('x',repmat(1:6,1,length(rejTrial)),'y',rejTrial2(:),'ymin',repmat(0,length(rejTrial2(:)),1),'ymax',repmat(100,length(rejTrial2(:)),1));
    g.stat_boxplot();
    g.set_names('x','condition','y','% of rejected epochs');
    if exp==1
        g.set_title('E1shortDC-LR');
    else
        g.set_title('E2longDC-LR');
    end
    g.draw(); 
    saveas(gcf,['figures' filesep 'E' num2str(exp) 'rejEpochsLR.pdf'])
   
    clear g rejTrial3;
    figure;
    rejTrial3 = rejTrial(:,7:12)';
    g = gramm('x',repmat(1:6,1,length(rejTrial)),'y',rejTrial3(:),'ymin',repmat(0,length(rejTrial3(:)),1),'ymax',repmat(100,length(rejTrial3(:)),1));
    g.stat_boxplot();
    g.set_names('x','condition','y','% of rejected epochs');
    if exp==1
        g.set_title('E1shortDC-SR');
    else
        g.set_title('E2longDC-SR');
    end
    g.draw(); 
    saveas(gcf,['figures' filesep 'E' num2str(exp) 'rejEpochsSR.pdf'])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create new sbjprediction with average ref for filteredWave and noiseWave
% ATTENTION: other fields have not been changed!
clearvars;
dataPath = {'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E1LRshortDC\',...
    'C:\Users\Marlene\Documents\JUSTIN\data\AMshortlong\E2LRlongDC\'};

for exp=1:length(dataPath)
   clear sbj;
   fileN = dir([dataPath{exp} 'sbjprediction.mat']);
   load([dataPath{exp} fileN.name]);
   for ss=1:length(sbj)
       for cond=1:size(sbj,2)
           averageAmp = (mean(sbj(ss,cond).data.filteredWave(:,:)));
           sbj(ss,cond).data.filteredWave = sbj(ss,cond).data.filteredWave - averageAmp;
           averageNoise = (mean(sbj(ss,cond).data.noiseWave(:,:)));
           sbj(ss,cond).data.noiseWave = sbj(ss,cond).data.noiseWave - averageNoise;           
       end
   end
   save([dataPath{exp} 'predAvRef'],'sbj');
end
