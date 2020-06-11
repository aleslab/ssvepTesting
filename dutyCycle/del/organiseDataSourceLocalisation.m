
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions

dirData = '/Users/marleneponcet/Documents/data/dutyCycle/';
cd(dirData)
listAxx = dir([dirData 'Axx' filesep '*.mat']);
newFolder = 'dataSourceLoc';

mkdir(newFolder); % create folder in which each sbj folder will be created
freqCond = {'F10_6','F5_3','F2_7'};

motion(1).cond = 18;
motion(2).cond = [17 19 22];
motion(3).cond = [16 20 21];
motion(1).label = 3;
motion(2).label = [2 4 3];
motion(3).label = [1 5 3];

for freq = 1:3
    % create separate folder for each freq condition
    mkdir([newFolder filesep freqCond{freq}])
    
    for ff=1:length(listAxx) 
        fprintf(['Processing S' num2str(ff) '\n'])
        
        % create sbj folder
        mkdir([newFolder filesep freqCond{freq} filesep listAxx(ff).name(end-6:end-4)])
        % create Exp folder containing Axx files
        axxFolder = [newFolder filesep freqCond{freq} filesep listAxx(ff).name(end-6:end-4) filesep 'Exp_Biosemi_128ch'];
        mkdir(axxFolder)
        % load Axx
        load ([dirData 'Axx' filesep listAxx(ff).name])
        
        ll=0;
        for cc=1+5*(freq-1):5*freq % do not include moving condition
            clear uniAxx;
            ll = ll+1;
            % double the wave for the single stim condition so that it fits
            % with the motion condition
            Axx(cc).wave = repmat(Axx(cc).wave,1,2);
            Axx(cc).nt = Axx(17+freq).nt; % nt from moving cond
            Axx(cc).time = Axx(17+freq).time; % time from moving cond
            % reformat Axx to be consistent with mrcDataViewer
            uniAxx = uniformAxx(Axx(cc));
            % save each condition separately
            save([axxFolder filesep 'Axx_c' num2str(ll,'%02d')],'-struct','uniAxx')
        end
        
        % add moving conditions
        for mm = 1:length(motion(freq).cond)
            clear uniAxx
            uniAxx = uniformAxx(Axx(motion(freq).cond(mm)));
            save([axxFolder filesep 'Axx_c' num2str(20+motion(freq).label(mm),'%02d')],'-struct','uniAxx')
        end
            
        
    end
    
end


%         % All the data that is loaded should be of the same dimensions
%         for cc=1:length(Axx)
%             sizeAxx(cc) = size(Axx(cc).wave,2);
%         end
%         if unique(sizeAxx) > 1
%             waveMax = max(sizeAxx);
%             for cc=1:length(Axx)
%                 if size(Axx(cc).wave,2) ~= waveMax
%                     nbCopy = waveMax / size(Axx(cc).wave,2);
%                     Axx(cc).wave = repmat(Axx(cc).wave,1,nbCopy);
%                     % need to replace nt as well
%                     Axx(cc).nt = waveMax;
%                     % I copy the time or should replace it?
%                     Axx(cc).time = repmat(Axx(cc).time,1,nbCopy);
%                 end
%             end
%         end