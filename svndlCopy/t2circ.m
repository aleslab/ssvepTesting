function [pVal, stdDev, confRadius, pT2, pChi] = t2circ(complexInput,alpha)
% t2circ - Calculate the T2circ statistic from Victor and Mast 1991
%function [pVal stdDev confRadius pT2 pChi] = tcirc(complexInput, [alpha])
%
%Input: compexInput: A vector or matrix of complex coefficients. In the
%                    case of a vector input returns scalar values. In the
%                    case of matrix input the calculation is done across
%                    the 1st dimension: Input d1 x d2 returns output 1 x d2.
%alpha:       [default=.05] Used for calculating confidence radius.
%
%Output:
%pVal:   Tcirc calculated pVal
%stdDev: Tcirc estimate of circular standard deviation. Pooled over
%        real/imag and assuming zero covariance.
%confRadius: 100*(1-ALPHA)% confidence radius
%pT2:    Hotelling T2 calculated pVal. This allows covariance. Has less
%        power than tCirc. *NOTE!* Not supported for matrix inputs. 
%
%(the following output are not to be used) 
%pChi: Don't use this!. This is a chi^2 approximate of the F table for the
%      T2. Implemented purely for internal diagnostics
%
%

%JMA
if isreal(complexInput)
	error('Vector Not Complex!')
end

if nargin<2
    alpha = .05;
end

if isvector(complexInput)
    vectorLength = length(complexInput);
else
    vectorLength = size(complexInput,1);
end;
    
M=vectorLength;

if (2*M-2)<=0
    error('Vector needs at least 2 elements');
end

realPart = real(complexInput);
imagPart = imag(complexInput);

realVar = var(realPart,0);
imagVar = var(imagPart,0);

%Equivalent to equation 1 in Victor and Mast
Vindiv =(realVar+imagVar)/2;

%Equation 2 of Victor and Mast
%length of mean vector squared
%!Imporant: Assuming test against hypothetical population mean of 0!
Vgroup = (M/2)*abs(mean(complexInput)).^2;

T2Circ = (Vgroup./Vindiv);

pVal = 1-fcdf(T2Circ,2,2*M-2);

stdDev = sqrt(Vindiv);

%We're going to do a silly thing here and cache outputs from finv.m
%This is a really-really slow function.  But in most use cases we get the
%same number over and over (e.g. getting t2circ for each fft frequency the
%alpha and M don't change.  
%I'm being a bit lazy by not making better vectorization code.  
%Also by throwing away cache if alpha value changes. 
%Also not checking that this cache doesn't get stupid big.  -JMA

% persistent fInvCache requestedAlpha;
% 
% if isempty(requestedAlpha);
%     requestedAlpha = alpha;
% end
% 
% if isempty(fInvCache) || alpha~=requestedAlpha
%     fInvCache = sparse(1,M);
% end
% 
% if size(fInvCache,2)>M || fInvCache(1,M) ~=0 %Check for cache hit and use it if exists
%     fInvVal = fInvCache(1,M);
% else
%     fInvVal = finv(1-alpha,2,2*M-2);
%     fInvCache(1,M) = fInvVal;
% end

fInvVal = finv(1-alpha,2,2*M-2);
%Equivalent to equation 5:
confRadiusSquared=2/M * fInvVal*Vindiv;
%To get the radius take square root.
confRadius = sqrt(confRadiusSquared);

%Should we return Hotelling values that don't assume equal variance.
%Not updated for matrix inputs yet - jma
if nargout>3
    
    if ~isvector(complexInput)
        error('Calculating Hotelling T2 values not supported for matrix sized inputs')
    end
    
 realMatrix = [real(complexInput) imag(complexInput)];
% 
[n,p]=size(realMatrix);
% 
m=mean(realMatrix,1); %Mean vector from data matrix X.
S=cov(realMatrix);  %Covariance matrix from data matrix X.
% S=eye(p)*Vindiv;
T2=n*(m)*inv(S)*(m)'; %Hotelling's T-Squared statistic.
F=(n-p)./((n-1)*p).*T2;

v1=p;  %Numerator degrees of freedom. 
v2=n-p;  %Denominator degrees of freedom.
%Probability that null Ho: is true. Test using F distribution


pT2=1-fcdf(F,v1,v2);  
% 



if nargout==5
    %Probability that null Ho: is true. Test using Chi^2 distribution
    v=p; %Degrees of freedom.
    pChi=1-chi2cdf(T2,v);
end

end
