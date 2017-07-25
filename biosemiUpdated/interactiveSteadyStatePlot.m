function [] = interactiveTopoSpecPlot(cfg, Axx)
%function [] = interactiveTopoPlot(cfg,Axx)
%
% helpful help here
%



% prepare the layout, this should be done only once
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg);

freq = Axx.freq;
time = 1000*Axx.time; %TODO: Make time units more explicit. 
Amp = Axx.Amp;
Wave = Axx.Wave;

if ~isfield(Axx,'tCircPval')
    sigFreqs = [];
else
    sigFreqs = Axx.tCircPval<=.01;
end

%Set default plot to show first:
iElec = 1;
iFr = 1;
iT = 1;


figH= figure('units','normalized','outerposition',[0 1 .6 .5]); %Render a new figure.

%Default Topography is spec. 
%topoAx = subplot(10,1,1:8);
topoAx = axes('Parent',figH,'Position',[0.05 .6 .3 .3]);
topoH = plotTopo(squeeze(Amp(:,iFr)),cfg.layout);
colormap(hot);

colorbarH = colorbar('peer',topoAx,'WestOutside')

colorbarH.Label.String = 'Microvolts';

cpos = colorbarH.Position;
cpos(3) = 0.5*cpos(3);
colorbarH.Position = cpos;

set(gcf,'KeyPressFcn',@keyInput)


axis off;
elecVerts = get(topoH,'Vertices');
hold on;
markH = plot(elecVerts(iElec,1),elecVerts(iElec,2),'ko','markersize',10,'linewidth',2);
elecNumH = text(-.05,1.35,num2str(iElec));


%Setup the frequency domain plot
%sigFreqs = (handles.data.i1F1:handles.data.i1F1:handles.data.nFr-1)+1;
specAx = axes('Parent',figH,'Position',[.4 .65 .5 .27]);
drawSpec();

%Setup the time domain plot
waveAx = axes('Parent',figH,'Position',[.4 .1 .5 .27]);
timeLine = [];
butterflyH = [];
selectedLineH = [];
drawWave();


%Setup complex phasor plot
phasorAx = axes('Parent',figH,'Position',[0.05 .1 .3 .3]);
drawPhase();

set(topoH,'ButtonDownFcn',@clickedTopo)
set(topoAx,'ButtonDownFcn',@clickedTopo)




%set(topoAx,'ButtonDownFcn',@specUpdate)

    function drawSpec()
        
        axes(specAx)
        [barH sigH] = pdSpecPlot(freq,Amp(iElec,:),sigFreqs(iElec,:));
        title(specAx,['Frequency: ' num2str(freq(iFr)) ' Hz'])
        
        %Set up the function to call when the plots are clicked on
        set(barH,'ButtonDownFcn',@clickedSpec)
        set(sigH,'ButtonDownFcn',@clickedSpec)
        set(specAx,'ButtonDownFcn',@clickedSpec)
        
    end


    function clickedSpec(varargin)
        
       tCurPoint = get(specAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance iFr] = min(abs(freq-tCurPoint(1,1)));
        
        if iFr>=1 && iFr<=length(freq)

            set(topoH,'facevertexCData',Amp(:,iFr))
            caxis(topoAx,[0 max(abs(Amp(:,iFr)))])
            title(specAx,['Frequency: ' num2str(freq(iFr)) ' Hz'])            
            colormap(hot);            
            drawPhase();

        else
            disp('clicked out of bounds')
        end
    end

    function clickedTopo(varargin)
        
        tCurPoint = get(topoAx,'CurrentPoint');
        
        dist = bsxfun(@minus,tCurPoint(1,:),elecVerts);
        
        dist = sqrt(sum(dist.^2,2));
        
        %Get the index to the nearest clicked electrode
        [distance iE] = min(dist);
        
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        if iE>=1 && iE<=size(elecVerts,1),
            
            iElec = iE;
            delete(markH);
            markH = plot(elecVerts(iElec,1),elecVerts(iElec,2),'ko','markersize',15,'linewidth',2);
            set(elecNumH,'String',num2str(iElec));
            title(topoAx,[num2str(iElec) ': ' cfg.layout.label{iElec}]);
            drawSpec();
            drawWave();
            drawPhase();

        end
        
    end

 function drawWave()

        %cla(waveAx)
     
        axes(waveAx)
        butterflyH = plot(waveAx,time,Wave','-','color',[.5 .5 .5]);
        hold on;       
        delete(selectedLineH);
        selectedLineH = plot(waveAx,time,Wave(iElec,:),'k','linewidth',2);  
        
        yLims = 1.1*[-max(abs(Wave(:))) max(abs(Wave(:)))];
        axis(waveAx,[0 time(end) yLims])
        
        [axLim] = axis(waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        delete(timeLine);
        timeLine = line([time(iT) time(iT)],[yLo yHi],'linewidth',2,'buttondownFcn',@clickedWave);
        
        ylabel('Amplitude (uV)');
        xlabel('Time (ms)');
        %Set up the function to call when the plots are clicked on
        set(waveAx,'ButtonDownFcn',@clickedWave)
        set(selectedLineH,'ButtonDownFcn',@clickedWave)
        set(butterflyH,'ButtonDownFcn',@clickedWave)
        
        
    end

    function clickedWave(varargin)
        
        tCurPoint = get(waveAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance iT] = min(abs(time-tCurPoint(1,1)));
        
        [axLim] = axis(waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        
        axes(waveAx);
        if iT>=1 && iT<=size(Amp,2)
            
            set(topoH,'facevertexCData',Wave(:,iT))
            %           caxis(topoAx,[-max(abs(data(iT,:))) max(abs(data(iT,:)))])
            
            colormap(jmaColors('arizona'));
            caxis(topoAx,[-max(abs(Wave(:))) max(abs(Wave(:)))]);            
            
            set(timeLine,'XData',[time(iT) time(iT)],'YData',[yLo yHi]);
            title(waveAx,['Time: ' num2str(time(iT),4) ' ms']);
            
        else
            disp('clicked out of bounds')
        end
    end


    function drawPhase()
        axes(phasorAx);
        delete(allchild(phasorAx)); %pdPhasePlot uses weird plotting functions so all it's objects need to be cleared to delete the scale. 
        pdPhasePlot( complex( Axx.Cos(iElec,iFr), Axx.Sin(iElec,iFr)),Axx.tCircStdErr(iElec,iFr));
        
    end

    function keyInput(src,evnt)
        
        switch(lower(evnt.Key))
            case 'leftarrow'
                iFr = max(iFr-1,1);
            case 'rightarrow'
                iFr= min(iFr+1,length(freq));
                
                
        end
        
        
          set(topoH,'facevertexCData',Amp(:,iFr))
          caxis(topoAx,[0 max(abs(Amp(:,iFr)))])
          title(specAx,['Frequency: ' num2str(freq(iFr)) ' Hz'])
            
        
    end
        
end
