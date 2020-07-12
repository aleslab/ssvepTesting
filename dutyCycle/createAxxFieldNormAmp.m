% create new field in Axx for normalised amplitudes

% normalise the frequency amplitudes by the amplitude of a square wave
% the idea here is that if a square wave gives a high amplitude for low DC
% and low for high DC (as we have in the data) then after normalisation we
% should bserve a flat pattern. Importantly, the normalisation also
% controls for the difference in amplitude due to the number of harmonics
% that are included in the signal (ie there are more for a 2.5 Hz than for
% a 10 Hz)

% in practice: get the sum of amplitudes for all the harmonics up to 50Hz
% for each condition separately for the recorded data and the square
% function. Then data./square for normalising

clearvars;
% load data
dataDir = '/Volumes/Amrutam/Marlene/JUSTIN/DutyCycle/data/Axx/';
listData = dir([dataDir '*.mat']);
keepSbj = [1:5 7:11 13:15 17:18 20]; 
% reject S7 and S13 S17, S20 same as S03
% S1=pilote so reject numS-1 (6 12 16 19)

% load the amplitudes for a square function
load('squareFFT.mat')

for ss = 1:length(keepSbj)
    clear Axx;
    load([dataDir listData(keepSbj(ss)).name]);
    
    % determine the indexes of all the harmonics < 50 Hz
    for cond=1:length(Axx)
        if cond < 16
            Axx(cond).harmIdx = determineFilterIndices( 'nf1low49', Axx(cond).freq, Axx(cond).i1f1);
        else
            Axx(cond).harmIdx = determineFilterIndices( 'nf1low49', Axx(cond).freq, Axx(cond).i1f1*2-1);
        end
    end
    
    % sum the harmonics and normalise
    for cond = 1:length(Axx)
        sumSq = sum(sqAllFFT(Axx(cond).harmIdx,cond));
        for ch=1:Axx(cond).nchan
            sumData = sum(Axx(cond).amp(ch,Axx(cond).harmIdx));
            Axx(cond).normFq(ch) = sumData / sumSq;
        end
    end
    
    save([dataDir listData(keepSbj(ss)).name],'Axx')

end