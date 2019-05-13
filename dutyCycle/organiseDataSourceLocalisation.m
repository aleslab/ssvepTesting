
dirAxx = dir(['/Users/marleneponcet/Documents/data/dutyCycle/Axx' filesep '*.mat']);
newFolder = 'dataSourceLoc';

mkdir(newFolder); % create folder in which each sbj folder will be created 

for ff=1:length(dirAxx)
    % create sbj folder
    mkdir([newFolder filesep dirAxx(ff).name(end-6:end-4)])
    % create Exp folder containing Axx files
    axxFolder = [newFolder filesep dirAxx(ff).name(end-6:end-4) filesep 'Exp_Biosemi_128ch'];
    mkdir(axxFolder)
    % load Axx
    load (['Axx' filesep dirAxx(ff).name])

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
    
    for cc=1:length(Axx)
        clear uniAxx;
        % reformat Axx to be consistent with mrcDataViewer
        uniAxx = uniformAxx(Axx(cc));
    	% save each condition separately
        save([axxFolder filesep 'Axx_c' num2str(cc,'%02d')],'-struct','uniAxx')       
    end

end


