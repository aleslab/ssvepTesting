clear all;
nBoot=10000;
nSamp = 50;
effectSize = 1;
pTcirc = zeros(nBoot,1);
stdDev = zeros(nBoot,1);
vecMean= zeros(nBoot,1);
pT2 = zeros(nBoot,1);
pChi = zeros(nBoot,1);
pZ = zeros(nBoot,1);
pZPooled = zeros(nBoot,1);
pTcircTrue = zeros(nBoot,1);
stdDevTrue = zeros(nBoot,1);
pT2True = zeros(nBoot,1);
pChiTrue = zeros(nBoot,1);

disp(['Running ' num2str(nBoot) ' bootstrap samples']);
for iBoot=1:nBoot

    %For null distribution
    complexVector = complex(randn(nSamp,1),randn(nSamp,1));
    vecMean(iBoot) = mean(complexVector);
    [pTcirc(iBoot) stdDev(iBoot) pT2(iBoot) pChi(iBoot) ] = tcirc(complexVector);
%    [pZ(iBoot) pZPooled(iBoot)] = zGilles(complex(randn,randn),complexVector);
    
    %For true positive distribution
%    [pZTrue(iBoot) pZPooledTrue(iBoot)] = zGilles(complex(effectSize+randn,randn),complexVector);
    complexVector = complex(randn(nSamp,1)+effectSize,randn(nSamp,1));
    [pTcircTrue(iBoot) stdDevTrue(iBoot) pT2True(iBoot) pChiTrue(iBoot) ] = tcirc(complexVector);
    
end


%check pval calibration:
disp('------------------------------------------------------------------------------');
disp('Checking T-CIRC NULL pVal calibration using empirical distribution from bootstrap:')
disp(['p-value exceeds .1 : ' num2str([sum(pTcirc<.1)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str([sum(pTcirc<.05)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str([sum(pTcirc<.01)/nBoot]) ' % of bootstrap samples'])
n = hist(pTcirc,round(nBoot/100));
expectedN = nBoot/round(nBoot/100);
disp(['p-value percent error: ' num2str(std(n)/expectedN)])


disp('------------------------------------------------------------------------------');
disp('Checking Hotelling T2 NULL pVal calibration using empirical distribution from bootstrap:')

disp(['p-value exceeds .1 : ' num2str([sum(pT2<.1)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str([sum(pT2<.05)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str([sum(pT2<.01)/nBoot]) ' % of bootstrap samples'])

n = hist(pT2,round(nBoot/100));
expectedN = nBoot/round(nBoot/100);
disp(['p-value percent error: ' num2str(std(n)/expectedN)])

disp('------------------------------------------------------------------------------');
disp('Checking Chi^2 approximation of NULL pVal calibration using empirical distribution from bootstrap:')

disp(['p-value exceeds .1 : ' num2str([sum(pChi<.1)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .05 : ' num2str([sum(pChi<.05)/nBoot]) ' % of bootstrap samples'])
disp(['p-value exceeds .01 : ' num2str([sum(pChi<.01)/nBoot]) ' % of bootstrap samples'])

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




