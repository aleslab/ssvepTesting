clear all;
nBoot=1000;
nSamp = 50;
effectSize = 1;
%Pre allocate all the different values
pTcirc = zeros(nBoot,1);
stdDev = zeros(nBoot,1);
confRadius = zeros(nBoot,1);
confRadiusTrue = zeros(nBoot,1);
vecMean= zeros(nBoot,1);
vecMeanTrue = zeros(nBoot,1);
pT2 = zeros(nBoot,1);
pChi = zeros(nBoot,1);
pZ = zeros(nBoot,1);
pZPooled = zeros(nBoot,1);
pTcircTrue = zeros(nBoot,1);
stdDevTrue = zeros(nBoot,1);
confRadiusTrue = zeros(nBoot,1);
pT2True = zeros(nBoot,1);
pChiTrue = zeros(nBoot,1);

%Example of theorietical calculation of confidence interval given the known
%population values in the simulation:
M = nSamp;
Vindiv = 4;
alpha = .05;
confRadiusTheory = sqrt(2/M * finv(1-alpha,2,2*M-2)*Vindiv);

disp(['Running ' num2str(nBoot) ' bootstrap samples']);
for iBoot=1:nBoot

    %For null distribution
    complexVector = complex(sqrt(Vindiv)*randn(nSamp,1),sqrt(Vindiv)*randn(nSamp,1));
    vecMean(iBoot) = mean(complexVector);
    [pTcirc(iBoot) stdDev(iBoot) confRadius(iBoot) pT2(iBoot) pChi(iBoot) ] = t2circ(complexVector);
%    [pZ(iBoot) pZPooled(iBoot)] = zGilles(complex(randn,randn),complexVector);
    
    %For true positive distribution
%    [pZTrue(iBoot) pZPooledTrue(iBoot)] = zGilles(complex(effectSize+randn,randn),complexVector);
    complexVector = complex(randn(nSamp,1)+effectSize,randn(nSamp,1));
    vecMeanTrue(iBoot) = mean(complexVector);
    [pTcircTrue(iBoot) stdDevTrue(iBoot) confRadiusTrue(iBoot) pT2True(iBoot) pChiTrue(iBoot) ] = t2circ(complexVector);
    
end



%check pval calibration:
disp('------------------------------------------------------------------------------');
disp('Checking T-CIRC NULL pVal calibration using empirical distribution from bootstrap:')
disp(['p-value exceeds .1 : ' num2str(100*sum(pTcirc<.1)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str(100*sum(pTcirc<.05)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str(100*sum(pTcirc<.01)/nBoot) ' % of bootstrap samples'])
n = hist(pTcirc,round(nBoot/100));
expectedN = nBoot/round(nBoot/100);
disp(['p-value percent error: ' num2str(std(n)/expectedN)])


disp('------------------------------------------------------------------------------');
disp('Checking T-CIRC 95% confidence radius using empirical distribution from bootstrap:')

%Does the observed mean plus confidence interval include 0? 
%To test see if the length of observed mean is more than a confidence
%radius away from 0. 
exceedConfInt = sum(abs(vecMean)>confRadius)/nBoot;

disp(['True value is outside 95% radius in: ' num2str(100*exceedConfInt) ' % of bootstrap samples'])



disp('------------------------------------------------------------------------------');
disp('Checking Hotelling T2 NULL pVal calibration using empirical distribution from bootstrap:')

disp(['p-value exceeds .1 : ' num2str(100*sum(pT2<.1)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str(100*sum(pT2<.05)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str(100*sum(pT2<.01)/nBoot) ' % of bootstrap samples'])

n = hist(pT2,round(nBoot/100));
expectedN = nBoot/round(nBoot/100);
disp(['p-value percent error: ' num2str(std(n)/expectedN)])

disp('------------------------------------------------------------------------------');
disp('Checking Chi^2 approximation of NULL pVal calibration using empirical distribution from bootstrap:')

disp(['p-value exceeds .1 : ' num2str(100*sum(pChi<.1)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str(100*sum(pChi<.05)/nBoot) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str(100*sum(pChi<.01)/nBoot) ' % of bootstrap samples'])

n = hist(pChi,round(nBoot/100));
expectedN = nBoot/round(nBoot/100);
disp(['p-value percent error: ' num2str(std(n)/expectedN)])


% disp('------------------------------------------------------------------------------');
% disp('Checking Z score NULL pVal calibration using empirical distribution from bootstrap:')
% 
% disp(['p-value exceeds .1 : ' num2str([sum(pZ<.1)/nBoot]) ' % of bootstrap samples'])
% disp(['p-value exceeds .05 : ' num2str([sum(pZ<.05)/nBoot]) ' % of bootstrap samples'])
% disp(['p-value exceeds .01 : ' num2str([sum(pZ<.01)/nBoot]) ' % of bootstrap samples'])
% 
% n = hist(pZ,round(nBoot/100));
% expectedN = nBoot/round(nBoot/100);
% disp(['p-value percent error: ' num2str(std(n)/expectedN)])
% 


% disp('------------------------------------------------------------------------------');
% disp('Checking Z score with pooled variance NULL pVal calibration using empirical distribution from bootstrap:')
% 
% disp(['p-value exceeds .1 : ' num2str([sum(pZPooled<.1)/nBoot]) ' % of bootstrap samples'])
% disp(['p-value exceeds .05 : ' num2str([sum(pZPooled<.05)/nBoot]) ' % of bootstrap samples'])
% disp(['p-value exceeds .01 : ' num2str([sum(pZPooled<.01)/nBoot]) ' % of bootstrap samples'])
% 
% n = hist(pZPooled,round(nBoot/100));
% expectedN = nBoot/round(nBoot/100);
% disp(['p-value percent error: ' num2str(std(n)/expectedN)])


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Power analysis')
disp('------------------------------------------------------------------------------');
disp('Checking T-CIRC Power using empirical distribution from bootstrap:')
disp(['Sample size: ' num2str(nSamp) ' Effect Size: ' num2str(effectSize)]);

disp(['hit rate at 0.1 : ' num2str(100*[sum(pTcircTrue<.1)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .05 : ' num2str(100*[sum(pTcircTrue<.05)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .01 : ' num2str(100*[sum(pTcircTrue<.01)/nBoot]) ' % of bootstrap samples'])

disp('------------------------------------------------------------------------------');
disp('Checking T-CIRC 95% confidence radius using empirical distribution from bootstrap:')


%Does the observed mean plus confidence interval include 0? 
%To test see if the length of observed mean is more than a confidence
%radius away from 0. 
%For true effect we subtract off the known population mean first. 
vecMeanTrueZeroed = vecMeanTrue - (effectSize+1i*0);
exceedConfIntTrue = 100*sum(abs(vecMeanTrueZeroed)>confRadiusTrue)/nBoot;
disp(['True value is outside 95% radius in: ' num2str(exceedConfIntTrue) ' % of bootstrap samples'])


disp('------------------------------------------------------------------------------');
disp('Checking Hotelling T2 Power using empirical distribution from bootstrap:')
disp(['Sample size: ' num2str(nSamp) ' Effect Size: ' num2str(effectSize)]);

disp(['hit rate at 0.1 : ' num2str(100*[sum(pT2True<.1)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .05 : ' num2str(100*[sum(pT2True<.05)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .01 : ' num2str(100*[sum(pT2True<.01)/nBoot]) ' % of bootstrap samples'])


disp('------------------------------------------------------------------------------');
disp('Checking Chi^2 Power using empirical distribution from bootstrap:')
disp(['Sample size: ' num2str(nSamp) ' Effect Size: ' num2str(effectSize)]);

disp(['hit rate at 0.1 : ' num2str(100*[sum(pChiTrue<.1)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .05 : ' num2str(100*[sum(pChiTrue<.05)/nBoot]) ' % of bootstrap samples'])
disp(['hit rate at .01 : ' num2str(100*[sum(pChiTrue<.01)/nBoot]) ' % of bootstrap samples'])



% disp('------------------------------------------------------------------------------');
% disp('Checking Z Power using empirical distribution from bootstrap:')
% disp(['Sample size: ' num2str(nSamp) ' Effect Size: ' num2str(effectSize)]);
% 
% disp(['hit rate at 0.1 : ' num2str(100*[sum(pZTrue<.1)/nBoot]) ' % of bootstrap samples'])
% disp(['hit rate at .05 : ' num2str(100*[sum(pZTrue<.05)/nBoot]) ' % of bootstrap samples'])
% disp(['hit rate at .01 : ' num2str(100*[sum(pZTrue<.01)/nBoot]) ' % of bootstrap samples'])
% 
% 
% disp('------------------------------------------------------------------------------');
% disp('Checking pooled variance Z Power using empirical distribution from bootstrap:')
% disp(['Sample size: ' num2str(nSamp) ' Effect Size: ' num2str(effectSize)]);
% 
% disp(['hit rate at 0.1 : ' num2str(100*[sum(pZPooledTrue<.1)/nBoot]) ' % of bootstrap samples'])
% disp(['hit rate at .05 : ' num2str(100*[sum(pZPooledTrue<.05)/nBoot]) ' % of bootstrap samples'])
% disp(['hit rate at .01 : ' num2str(100*[sum(pZPooledTrue<.01)/nBoot]) ' % of bootstrap samples'])
% 
% 




