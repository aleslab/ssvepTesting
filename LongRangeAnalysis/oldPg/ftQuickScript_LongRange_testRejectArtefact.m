% still have to figure out how to remove noise and bad channels
% A23 = Oz

clear cfg data;

addpath /Users/marleneponcet/Documents/Git/fieldtrip
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
ft_defaults
cfg.dataset   =  '/Users/marleneponcet/Documents/LongRangeSSVEP/LongRangeS01.bdf';

% cfg.dataset   =  'D:\StAndrews\LongRangeSSVEP\LongRangeS01.bdf'
% addpath D:\GitHubRepo\fieldtrip
% addpath D:\GitHubRepo\ssvepTesting\svndlCopy
% addpath D:\GitHubRepo\ssvepTesting\biosemiUpdated

% cfg.dataset   = '/Volumes/Amrutam/Marlene/LongRangeSSVEP/LongRangeS01.bdf'
% addpath /Volumes/Amrutam/Marlene/Git/fieldtrip
% addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/svndlCopy
% addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated

% define trials then preprocessing or the other way around?? Faster this
% way.. 

% read the behavioural file
load('longRange_01__20170727_153334.mat')
tableData = struct2table(experimentData);


cfg.trialdef.preStimDuration = experimentData(1).trialData.preStimDuration;
cfg.trialdef.nbTotalCycles = experimentData(1).trialData.nbTotalCycles;
cfg.i1f1 = experimentData(1).stimTagFreq;

cfg.channel = 1:128;

cfg.trialdef.bitmask = 2^9-1; %Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.layout = 'biosemi128.lay';
cfg.trialfun = 'lock2ssvep_LongRange_V2'; 
[cfg] = ft_definetrial(cfg);

cfg.preproc.lpfilter  = 'no';
cfg.preproc.lpfreq        = 80;
cfg.preproc.demean        ='yes';
cfg.preproc.reref         = 'no'; 
cfg.preproc.refchannel    = 'A1'; % A3 = CPz / use Cz = A1
[data] = ft_preprocessing(cfg); 

% data = ft_redefinetrial(cfg,data);



%%%Next resample the data so we have an integer number of samples per cycle.
cfg.newFs = 85*8; %Integer number of samples per monitor refresh
data = resample_steadystate_test(cfg,data);
% because of reseampling, should change data.trialinfo
data.trialinfo(:,3) = size(data.time{1,1},2) / 5; % 5 cycles / epoch

%DO EPOCH ARTIFACT REJECTION HERE
cfg.layout = 'biosemi128.lay';
cfg.method = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg, data);

cfg.method = 'summary';
cfg.keepchannel = 'repair'; % had to modify ft_rejectvisual line 343 so that layout was taken into account
[cleanData] = ft_rejectvisual_modif(cfg, data);
% 340 trials
% 1, 3, 5, 6, 10, 13, 14, 19, 20, 28, 36, 43, 48, 52, 59, 64, 68, 73, 78, 80, 82, 83, 85, 86, 89, 93, 103, 106, 109, 110, 112, 116, 117, 122, 123, 124, 125, 126, 127, 128, 130, 135, 137, 142, 144, 145, 147, 149, 150, 152, 153, 154, 155, 157, 158, 162, 163, 164, 165, 166, 168, 169, 173, 175, 176, 183, 187, 189, 190, 191, 192, 194, 195, 196, 197, 199, 202, 203, 207, 213, 220, 222, 223, 226, 227, 229, 231, 242, 245, 247, 252, 253, 254, 256, 258, 260, 262, 270, 276, 281, 282, 284, 286, 287, 289, 292, 295, 305, 309, 313, 315, 316, 319, 320, 321, 322, 323, 324, 327, 330, 331, 334, 337, 338, 342, 349, 351, 352, 355, 356, 361, 363, 364, 367, 369, 370, 373, 379, 382, 383, 388, 389, 394, 395, 396, 397, 400, 401, 405, 406, 407, 409, 412, 413, 414, 418, 420, 421, 423, 425, 427, 428, 430, 431, 433, 438, 440, 442, 444, 445, 446, 448, 452, 455, 456, 457, 458, 459, 460, 461, 462, 464, 465, 467, 470, 472, 474, 475, 476, 478, 479, 481, 483, 485, 486, 489, 499, 502, 505, 506, 509, 514, 517, 518, 524, 527, 529, 530, 531, 533, 535, 536, 538, 540, 541, 542, 543, 544, 550, 552, 554, 557, 558, 559, 560, 561, 562, 565, 566, 569, 577, 578, 580, 582, 584, 585, 586, 587, 592, 598, 604, 607, 608, 610, 612, 615, 616, 619, 622, 623, 628, 632, 636, 639, 642, 645, 650, 651, 655, 657, 659, 663, 668, 672, 674, 676, 677, 678, 679, 681, 682, 683, 684, 687, 690, 692, 693, 694, 699, 700, 701, 702, 703, 707, 712, 713, 714, 717, 722, 724, 729, 731, 733, 734, 739, 743, 749, 751, 755, 759, 763, 772, 774, 779, 786, 789, 793, 795, 798, 803, 807, 809, 811, 820, 821, 826, 831, 834, 835, 841, 843, 845, 848, 851, 852, 855, 856, 857, 858, 861, 865, 871, 878, 881, 883, 889, 890, 894, 897, 900

cfg.artfctdef.threshold.range = 150;
cfg.artfctdef.threshold.bpfilter  = 'no';
cfg.continuous = 'no';
cfg.feedback = 'yes';
[cfg, artifact] = ft_artifact_threshold(cfg); % 311 artefacts

cfg = ft_databrowser(cfg,cleanData); % it shows flat for the reference channel (A1) but there is data in there...
cfg.artfctdef.reject  = 'complete';
cleandata = ft_rejectartifact(cfg,cleanData); % +29
% 531/585 trials kept
save('cleandata_A13C7', 'cleandata')

cfg = [];
cfg.vartrllength = 2;
[timelockEpoch] = ft_timelockanalysis(cfg, cleandata);

%remove the extra channels. 
cfg.channel =  {'all','-EXG1', '-EXG2', '-EXG3','-EXG4','-EXG5','-EXG6','-EXG7','-EXG8', '-Status'};
dataEeg = ft_selectdata(cfg,cleandata);
%Do the steadystate analysis
cfg.trimdur = 0;
[Axx] = ft_steadystateanalysis(cfg, cleandata);


% %Now do a quick spectrum of the 12th channel:
% figure;
% pdSpecPlot(Axx.freq(2:80),Axx.Amp(12,2:80)',Axx.tCircPval(12,2:80)<.05)
% title('Channel 12')
% 
% %Now setup an interactive plot
% figure;clf;
% interactiveTopoSpecPlot(cfg,Axx)
% interactiveSteadyStatePlot(cfg,Axx)


% %%% Power Spectrum
% fftDat = fft(timelockEpoch.avg(1:128,:)'); % for the 32 channels
% % spectrumDat = abs(2*fftDat(2:floor(end/2))/length(fftDat));
% spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric). 
% % Also transforms the matrix in a single vector. To keep the channels
% % separate (2D matrix): spectrumDat = abs(fftDat(2:floor(end/2),:));
% spectrumDat = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
% freqs = [1:length(spectrumDat)] * 1/(timelockEpoch.time(end)-timelockEpoch.time(1));
% figure; hold on; plot(freqs(freqs<80),spectrumDat(freqs<80),'r');
% pdSpecPlot(freqs(1:80),spectrumDat(1:80),[]);
% 
% spectrumChan = abs(fftDat(2:floor(end/2),16)); % for the specified channel, here 16=Oz?
% spectrumChan = spectrumChan*2/length(spectrumChan);
% figure; hold on; plot(freqs(freqs<80),spectrumChan(freqs<80),'r');
% pdSpecPlot(freqs(1:80),spectrumChan(1:80),[]);


% %%%%%%%%%%%%%%%%%%
% %%%% depending on condition (saved in data.trialinfo)
% %%%%%%%%%%%%%%%%%%
% 
% allcond = unique(data.trialinfo(:,1));
% for cond = 1 : length(allcond)
%     cfg = [];
%     cfg.trials = find(data.trialinfo(:,1) == allcond(cond));
%     condAverage(cond).data = ft_timelockanalysis(cfg, data);
%     fftDat = fft(condAverage(cond).data.avg(23,:)'); 
%     spectrumDat = abs(fftDat(2:floor(end/2))); % discard the second half of the signal (since it's symetric).
%     condSpect(cond,:) = spectrumDat*2/length(fftDat); % estimate of the power amplitude = 2*A/N because N sample data amplitude A results in a DFT peak of AN/2
% end
% 
% freqs = [1:length(spectrumDat)] * 1/(condAverage(cond).data.time(end)-condAverage(cond).data.time(1));
% 
% figure; hold on; 
% for cond = 1 : length(allcond)
%     subplot(3,3,cond)
%     plot(freqs(freqs<40),condSpect(cond,freqs<40),'r');
% end



%%%%% follow up
cfg.layout = 'biosemi128.lay';
cfg.channel = {'all'};
allcond = unique(cleandata.trialinfo(:,1));
for cond=1:length(allcond)
    cfg.trials = find(cleandata.trialinfo(:,1) == allcond(cond));
%     cfg.channel =  {'all','-GSR1', '-GSR2', '-Erg1','-Erg2','-Resp','-Plet','-Temp', '-Status'};
    % dataEeg = ft_selectdata(cfg,data);
    [Axx(cond)] = ft_steadystateanalysis(cfg, cleandata);
%     cfg.channel = {'all'};
end
cfg = [];
cfg.layout = 'biosemi128.lay';
save('AxxS01', 'Axx','cfg')
interactiveSteadyStatePlot(cfg,Axx)
