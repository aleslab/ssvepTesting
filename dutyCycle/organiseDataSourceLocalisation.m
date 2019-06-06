
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/commonFunctions

dirData = '/Users/marleneponcet/Documents/data/dutyCycle/';
cd(dirData)
listAxx = dir([dirData 'Axx' filesep '*.mat']);
newFolder = 'dataSourceLoc';

mkdir(newFolder); % create folder in which each sbj folder will be created
freqCond = {'F10_6','F5_3','F2_7'};

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
        
        % copy wave if it is not the same size as the other conditions
        for cc=1:length(Axx)
            sizeAxx(cc) = size(Axx(cc).wave,2);
        end
        if unique(sizeAxx) > 1
            waveMax = max(sizeAxx);
            for cc=1:length(Axx)
                if size(Axx(cc).wave,2) ~= waveMax
                    nbCopy = waveMax / size(Axx(cc).wave,2);
                    Axx(cc).wave = repmat(Axx(cc).wave,1,nbCopy);
                    % need to replace nt as well
                    Axx(cc).nt = waveMax;
                    % I copy the time or should replace it?
                    Axx(cc).time = repmat(Axx(cc).time,1,nbCopy);
                end
            end
        end
        
        for cc=1+5*(freq-1):5*freq % do not include moving condition
            clear uniAxx;
            % reformat Axx to be consistent with mrcDataViewer
            uniAxx = uniformAxx(Axx(cc));
            % save each condition separately
            save([axxFolder filesep 'Axx_c' num2str(cc,'%02d')],'-struct','uniAxx')
        end
        
    end
    
end
