clear cfg;
% To facilitate data-handling and distributed computing you can use
cfg.dataset   =  '/Users/ales/data/ssvepTesting/SSVEPtest3.bdf'

%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.padding      = length (in seconds) to which the trials are padded for filtering (default = 0)
%   cfg.padtype      = string, type of padding (default: 'data' padding or
%                      'mirror', depending on feasibility)
%   cfg.continuous   = 'yes' or 'no' whether the file contains continuous data
%                      (default is determined automatic)
%

cfg.trialdef.bitmask = 2^9-1;%Values to keep.
cfg.trialdef.condRange = [101 165];
cfg.trialdef.ssvepTagVal = 1;
cfg.trialdef.epochLength = 2; 
cfg.trialfun = 'lock2SsvepTag'; 

[cfg] = ft_definetrial(cfg)

cfg.demean        ='yes';
cfg.reref         = 'yes';
cfg.refchannel    = {'A1'};
[data] = ft_preprocessing(cfg)

[timelock] = ft_timelockanalysis(cfg, data)
 

Axx.Wave = timelock.avg(1:32,:)';
Axx.nT = size(Axx.Wave,1);
Axx.nFr = round(size(Axx.Wave,1)/2);
dft = dftmtx(Axx.nT);

dftDat = dft*Axx.Wave;
dftDat = dftDat(1:Axx.nFr,:);

Axx.dFHz = (data.hdr.Fs)/Axx.nT;

% Stuff tends to be multiples of 1 Hz. But we don't have that info
% here so we are just going to set the index to 1 to make nF1 be all
% freqs
Axx.i1F1       = 1;
Axx.i1F2       = 0;

Axx.Amp = abs(2*(dftDat/Axx.nFr));
Axx.Cos = 2*real(dftDat)/Axx.nFr;
Axx.Sin = -2*imag(dftDat)/Axx.nFr;



freqs = 0:Axx.dFHz:(Axx.dFHz*(Axx.nFr-1));

pdSpecPlot(freqs(1:80),Axx.Amp(1:80,12)',[])

