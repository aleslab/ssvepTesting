%%% simulation of EEG response across harmonics for different duty cycles
clearvars; close all;
ff=2.5; % fondamental frequency, min is 1.25 use multiple of 1.25


%% Generate onset/offset response (i.e. odd and even harmonics)
nT=512; % has to be divisible by 64
wtSize = 1.2; % size of time window = 1.2 sec
t=linspace(0,wtSize,nT); 

wList= [1:1:12]'; % frequencies used for the ERP
nF = length(wList);
wAmp = [linspace(.05,1,nF/2) linspace(1,.25,nF/2) ]'; % define the amplitude response change with frequencies (peak around 8Hz)
wAmp(2:2:end) = wAmp(2:2:end)*1;
%%%%%%%%%%%%% what is the phase for? useful for waveform, less for
%%%%%%%%%%%%% spectrum, looks better with the multiplication
wPhase = (-1*wList)+repmat([3*pi/2 4*pi/3]',nF/2,1); % freq + phase of odd and even harmonics [sinewave = sin(t+theta)]
% wPhase = (-1*wList)+repmat([0*pi/2 0*pi/2]',nF/2,1);

wAmp = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;


%% Make a fourier basis set for an easy waveform creator:
allDC = [0.125 0.25 0.5 0.75 0.875];

for dc=1:length(allDC)
    dutyCycle= allDC(dc);
    clear waveForm

    fourierBasis = dftmtx(nT);
    invFourier  = conj(fourierBasis)/nT;
    waveForm = wAmp'*fourierBasis(wList+1,:);
    waveForm = real(waveForm);
    
    % 1.25Hz
    if ff == 1.25
        waveForm = waveForm + circshift(waveForm,round(nT*dutyCycle)); % On + Off responses
    else
        % 2.5 Hz
        if ff==2.5
            nbCy = 2;
            %     figure;plot(circshift(circshift(waveForm,nT/2),round(nT/2*dutyCycle)))
            %     figure;plot(waveForm + circshift(waveForm,round(nT/2*dutyCycle)) + circshift(waveForm,nT/2) + circshift(circshift(waveForm,nT/2),round(nT/2*dutyCycle)));
            %         waveForm = waveForm + circshift(waveForm,round(nT/2*dutyCycle)) + circshift(waveForm,nT/2) + circshift(circshift(waveForm,nT/2),round(nT/2*dutyCycle));
        elseif ff==5
            nbCy = 4;
        elseif ff==10
            nbCy = 8;
        end
        % ON responses
        for mm=1:nbCy-1
            waveForm = waveForm + circshift(waveForm,mm*nT/nbCy);
        end
%         figure;plot(waveForm)
%         figure;plot(circshift(waveForm,(nT/nbCy)*dutyCycle))
        % OFF responses
        waveForm = waveForm + circshift(waveForm/4,(nT/nbCy)*dutyCycle);
    end
    
    waveForm = waveForm./max(abs(waveForm));
    
    wAmp2= nF*(abs(invFourier(wList+1,:)*waveForm'));
    
    %Make spec plot
    freqs = ff*(0:.5:20); % spectrum for fq until 20Hz
    amps = 0*ones(size(freqs));
    amps(2*wList+1)=abs(wAmp2);
    
    
    %% figures
    %timecourse panel
    figure(100)
    subplot(2,5,dc)
    plot(t,waveForm,'k');
    xlabel('Time (seconds)')
    ylabel('Amplitude');
    ylim([-1.2 1.2])
    xlim([0 wtSize])
    title([num2str(ff) 'Hz ' num2str(dutyCycle*100) '%'])

    subplot(2,5,dc+5)
    if exist('pdSpecPlot')
        pdSpecPlot(freqs,amps,[]);
    else
        bar(freqs,amps);
    end
    xlim([0 30])
    xlabel('Frequency (Hertz)')
    ylabel('Amplitude');
    title([num2str(ff) 'Hz ' num2str(dutyCycle*100) '%'])
end
set(100, 'Position', [0 0 1000 500])
% saveas(100,['figures' filesep 'simulHarmonics' num2str(ff) 'Hz'],'png')
