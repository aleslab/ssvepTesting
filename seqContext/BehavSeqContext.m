% analyse behavioural results for seqContext
% detection of white dots (0-3) on the black vertical bars
% filesDir = 'C:\Users\Marlene\Documents\dataStAndrews\seqContext\behaviour\';
filesDir = '/Users/marleneponcet/Documents/data/seqContext/originalData/';
allfiles = dir([filesDir '*.mat']);

for ff = 1:length(allfiles)
    load([filesDir allfiles(ff).name])
    
    % remove any invalid trial from the data
    toRemove = [];
    for ll=1:length(experimentData)
        if experimentData(ll).validTrial == 0
            toRemove = [toRemove; ll];
        end
    end
    keepIndexes = setdiff(1:length(experimentData), toRemove);
    tableData = struct2table(experimentData(keepIndexes));
    
    %%% analysis per condition 
    allCond = unique(tableData.condNumber);
    for cc=1:length(allCond)
        indexCond = find(tableData.condNumber == allCond(cc));
        corResp = [];
        for ic = 1: length(indexCond)
            if tableData.trialData(indexCond(ic)).nbDots == tableData.trialData(indexCond(ic)).response
                corResp(ic) = 1;
            else
                corResp(ic) = 0;
            end
        end
        acc(ff,cc) = mean(corResp);
    end
    
end

boxplot(acc)