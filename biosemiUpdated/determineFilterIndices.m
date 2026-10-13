function freqIndices = determineFilterIndices(filterName, freqValues, i1f1, i1f2)
% Select nF1/nF2 harmonics using one-based fundamental DFT-bin indices.
% Supports Odd, Even, Clean and LowNN suffixes, e.g. nF2CleanLow49.
% Clean excludes harmonics shared with the other tag. Supply i1f2 for
% nF2 filters and for removing f2 overlaps from nF1Clean. DC is excluded.
if nargin < 4, i1f2 = []; end
name = lower(filterName);
cutoff = Inf;
lowStart = regexp(name,'low\d+$','start','once');
if ~isempty(lowStart)
    cutoff = str2double(name(lowStart+3:end));
    name = name(1:lowStart-1);
end
if ~any(strcmp(name,{'nf1','nf1odd','nf1even','nf1clean','nf2','nf2odd','nf2even','nf2clean'}))
    error('determineFilterIndices:UnknownFilter','Unknown filter: %s',filterName);
end
fundamentalIdx = i1f1;
if strncmp(name,'nf2',3), fundamentalIdx = i1f2; end
step = fundamentalIdx-1;
harmonics = 1:floor((numel(freqValues)-1)/step);
if strcmp(name(4:end),'odd')
    harmonics = harmonics(mod(harmonics,2)==1);
elseif strcmp(name(4:end),'even')
    harmonics = harmonics(mod(harmonics,2)==0);
end
freqIndices = 1+harmonics*step;
if strcmp(name(4:end),'clean')
    otherIdx = i1f2;
    if strncmp(name,'nf2',3), otherIdx = i1f1; end
    if ~isempty(otherIdx)
        freqIndices = freqIndices(mod(freqIndices-1,otherIdx-1) ~= 0);
    end
end
freqIndices = freqIndices(freqValues(freqIndices) <= cutoff);
end
