function [ freqIndices ] = determineFilterIndices( filterName, freqValues, i1f1 )
%determineFilterIndices Finds the indices for a selected steady state filter
% [ freqIndices ] = determineFilterIndices( filterName, freqValues, i1f1 )
%
% filterName - A string specifying the filter:
%
%'nF1'     - ALl multiples of the steady state frequency
%'nF1Odd'  - All odd harmonics
%'nF1Even' - All even harmonics
%
% The above can be extended to include a lowpass cuttoff by adding:
% 'LowDD' where DD is an integer representing the highest frequency to
% include in the filter.
%
% freqValues - A vector including all the frequency values in the response.
%              E.g. [0 .5 1 1.5 2 2.5 ... 120.5];
%
% i1f1       - Index of the frequency of the fundamental in freqValues. 
%
% Example: 
%
% freqValues = 0:.5:20;
% i1f1 = 6; %Equal to 2.5 Hz from the freqValues(6);
% filterName = 'nF1Low10'; %Request all harmonics 10 or lower.
%
% freqIndices = determineFilterIndices(filterName,freqValues,i1f1);




%Filter definitions:
filterDef  = { 'nF1' 'ALl multiples of the steady state frequency';...
               'nF1Odd' 'All odd harmonics';...
               'nF1Even' 'All even harmonics';...
               'xxxLowFF', 'Same as others but only up to FF Hz'};


%use a regex to find if the filter spec requests lowpas.
if any(regexpi(filterName,'low\d*$')) 
    
    cutoffVal = regexpi(filterName,'\d*$','match');
    cutoffVal = str2double(cutoffVal);
    freqValues = freqValues(freqValues<=cutoffVal);
           
    
    lowpassChar = regexpi(filterName,'low\d*$');
    filterName = filterName(1:(lowpassChar-1)); %Now cutoff the lowpass part of the name
end

switch lower(filterName)
    
    %Note: the i1f1-1 is there to account for the DC component in the DFT
    case 'nf1'
        freqIndices = i1f1:(i1f1-1):length(freqValues);
    case 'nf1odd'
        freqIndices = i1f1:2*(i1f1-1):length(freqValues);
    case 'nf1even'
        freqIndices = (2*(i1f1-1):2*(i1f1-1):length(freqValues))+1;
    otherwise
        freqIndices = 2:length(freqValues); %Do not include DC
        warning('Unrecognized filter name');
        

end

