function [newData] = rerefData(data,ref)
% re-reference the data according to new ref channel
% data is Axx (elec x fq)
% ref is reference channel number

matVersion = version('-release');
matVersion = str2num(matVersion(1:4));
if matVersion == 2016
    if strcmp(matVersion(5),'a')
        matVersion = 2015;
    else
        matVersion = 2017;
    end
end

if matVersion<2016
    error('old matlab version not coded at the moment')
    % DataCzRef = bsxfun(@minus, data, CzRef)
else
    newData = data;
    newData.cos = data.cos-data.cos(ref,:);
    newData.sin = data.sin-data.sin(ref,:);
    newData.wave = data.wave-data.wave(ref,:);
    newData.amp = sqrt(newData.cos.*newData.cos+newData.sin.*newData.sin);
end

end