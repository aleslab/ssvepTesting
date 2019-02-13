function [ filteredWave ] = filterSteadyState( cfg, steadyState )
%filterSteadyState Reconstructs filtered waveform from chosen frequencies
% [ filteredWave ] = filterSteadyState( cfg, steadyState )   
%
%
%
%activeFreq is a matrix same size as Amp/Sin/Cos selecting active elements
%Create the component matrix
%First determine the length of the fourier transform components.

if ~isfield(steadyState, 'ndft') %If we're not told the length, have to try and recreate it. 
    
    ndft = 2*(steadyState.nfr-1); %Close but may be off by 1
    
    if mod(ndft,steadyState.nt)~=0,
        ndft = ndft +1;
        if mod(ndft,steadyState.nt)~=0,
            error('Failed to guess nDft. Are nT and nFr compatible?');
        end
    end
else
    ndft = steadyState.ndft;
end

%Now create the component matrix:
%Note: The help for dftmtx says the inverse transform is the complex
%conjugate divided by n.  This is true for the matrix multiplication.
%However we don't story the data that way. We store the the data coefficients in a
%way that makes sense for interpreting eeg data.  That is already
%normalized by nDft and the sin coefficient, which is the imaginary part
%has already been sign inverted.  Therefore we need to make our own inverse
%dftmtx. 
dft = dftmtx(ndft);
dft = dft(:,1:steadyState.nfr); %Just grab the frequency components we need. Do not include DC by default. 
filteredWave = NaN(size(steadyState.wave));


for iChan = 1:steadyState.nchan,
    
    selFr = cfg.activeFreq(iChan,:);
    %Recon: Wave = real(theta)*Cos(theta) - imag(theta)*Sin(theta)
    %We're going to recon the full fourier transform sized waveform.
    %Then reduce it down the the single cycle length. This way we can
    %include the noise bands in the waveform if we want, or we can just use
    %the cycle components. 
    waveRecon = real(dft(:,selFr))*steadyState.cos(iChan,selFr)' - imag(dft(:,selFr))*steadyState.sin(iChan,selFr)';    
    waveRecon=reshape(waveRecon,steadyState.nt,[])';
    filteredWave(iChan,:)=mean(waveRecon,1);
end



end

