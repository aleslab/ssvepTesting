function newAxx = uniformAxx(Axx)
% make field names uniform in Axx from fieldtrip to mrcDataViewer 
% and transpose matrix elec x fq to fq x elec
% there might be other discrepancies not corrected here... 
% also it does not check for differences in nb of fields!


% Modify field names
oldFields = fieldnames(Axx);
newFields = { 'nCh','fsample','nFr','ndft','nT','dTms','time','Wave',...
    'Amp','Sin','Cos','freq','dFHz','i1F1','pval','stderradius',...
    'confradius','label','dimord','wavedimord','elec','cfg','condLabel'};

for k=1:numel(oldFields)
    [newAxx.(newFields{k})] = deal(Axx.(oldFields{k})) ;
end

% transpose matrices 
newAxx.Wave = newAxx.Wave';
newAxx.Amp = newAxx.Amp';
newAxx.Sin = newAxx.Sin';
newAxx.Cos = newAxx.Cos';
newAxx.pval = newAxx.pval';
newAxx.stderradius = newAxx.stderradius';
newAxx.confradius = newAxx.confradius';

% add field i1F2
% no second harmonic but required
newAxx.i1F2 = 0;
% first harmonic should be -1 since DC already taken care of in the pg
newAxx.i1F1 = newAxx.i1F1-1;
