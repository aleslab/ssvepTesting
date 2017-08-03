function [] = interactiveSteadyStatePlot(cfg, steadyState)
%function [] = interactiveTopoPlot(cfg,steadyState)
%
% helpful help here
%


%Setup default condition choices.
configOptions.idxCondA = 1;
configOptions.idxCondB = 1;
if length(steadyState)>=2, %If 2 or more conditions default to showing cond 2. 
    configOptions.idxCondB = 2;
end

configOptions.pValThresh = Inf;

%TODO: Clarify and check this code. 
% prepare the layout, this should be done only once
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg);



%
%plotConfig - A structure containing the different options to use for
%             displaying data


%Using condData to hold th 2 conditions to plot. 
condData(1) = initPlotData(steadyState(configOptions.idxCondA));
condData(2) = initPlotData(steadyState(configOptions.idxCondB));
condData(1).selectedCond = configOptions.idxCondA;
condData(2).selectedCond = configOptions.idxCondB;


condData(1).color = [1 0 0];
condData(2).color = [0 0 1];

iElec = 1;

%Initialize figure;
figH= figure('units','normalized','outerposition',[0 1 .8 1]); %Render a new figure.




%%%%% Setup plot option GUI controls
% Create pop-up menu

%Quickie Just label conditions with numbers:
conditionList = cellstr(num2str([1:length(steadyState)]'));

pValList = {'None','0.05','0.01','1e-3','1e-6'};
filterList = {'None','nF1','nF1Odd','nF1Even','nF1Low49'};
    
uicontrol('Style', 'text','units','normalized',...
    'String','A',...
    'Position', [.9 .975 .05 .02])

uicontrol('Style', 'text','units','normalized',...
    'String','B',...
    'Position', [.95 .975 .05 .02])

uicontrol('Style', 'text','units','normalized',...
    'String','Condition: ',...
    'Position', [.85 .94 .05 .02])

selCondAPopupH = uicontrol('Style', 'popup','units','normalized',...
    'String',conditionList,'Value',configOptions.idxCondA,...
    'Position', [.9 .95 .05 .01],...
    'Callback', {@selCond,1});

selCondBPopupH = uicontrol('Style', 'popup','units','normalized',...
    'String',conditionList,'Value',configOptions.idxCondB,...
    'Position', [.95 .95 .05 .01],...
    'Callback', {@selCond,2});


uicontrol('Style', 'text','units','normalized',...
    'String','p Threshold:',...
    'Position', [.82 .9 .08 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',pValList,'Value',1,...
    'Position', [.9 .9 .1 .02],...
    'Callback', @selPVal);

uicontrol('Style', 'text','units','normalized',...
    'String','Time Filter:',...
    'Position', [.81 .85 .06 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',filterList,'Value',1,...
    'Position', [.87 .85 .07 .02],'fontsize',8,...
    'Callback', {@selFilter,1});

uicontrol('Style', 'popup','units','normalized',...
    'String',filterList,'Value',1,...
    'Position', [.93 .85 .07 .02],'fontsize',8,...
    'Callback', {@selFilter,2});



infoPanelH = uipanel('Title','Info','FontSize',12,...
             'BackgroundColor','white',...
             'Position',[.82 .4 .16 .4],'visible','off');


%%%%%%% Setup data plot axes.

%Default Topography is spec. 
%topoAx = subplot(10,1,1:8);
condData(1).topoAx = axes('Parent',figH,'Position',[0.05 .6 .3 .3]);
initTopo(1);

condData(2).topoAx = axes('Parent',figH,'Position',[0.05 .25 .3 .3]);
initTopo(2);

%Setup the frequency domain plot
%For condition A
condData(1).specAx = axes('Parent',figH,'Position',[.4 .82 .4 .12]);
drawSpec(1);

%For condition B
condData(2).specAx = axes('Parent',figH,'Position',[.4 .57 .4 .12]);
drawSpec(2);

%Setup the time domain plot
condData(1).waveAx = axes('Parent',figH,'Position',[.4 .35 .4 .12]);
drawWave(1);

condData(2).waveAx = axes('Parent',figH,'Position',[.4 .1 .4 .12]);
drawWave(2);



%Setup complex phasor plot
phasorAx = axes('Parent',figH,'Position',[0.85 .1 .14 .14]);
drawPhase();


%%%% Setup info pane information
drawInfoPane()


%set(topoAx,'ButtonDownFcn',@specUpdate)


    function plotData = initPlotData(steadyState)
    %This function does the setup needed to take the steadyState fields and
    %make them into values needed for plotting.
        
    plotData = steadyState;
    
    %!!!!!! This is a really, really stupid line.  Just relying on stupidy
    %being obvious on the plot if units are off by x1000
    plotData.Time = 1000*steadyState.time; %TODO: Make time units more explicit.        
    
    if ~isfield(steadyState,'tCircPval')
        plotData.sigFreqs = [];
    else
        plotData.sigFreqs = steadyState.tCircPval<=configOptions.pValThresh;
    end
    
    %Set default plot options        
    plotData.iElec = 1;
    plotData.iFr = 1;
    plotData.iT = 1;
    
    %Handles for plot elements
    plotData.timeLine = [];
    plotData.butterflyH = [];
    plotData.selectedLineH = [];
        
    end

    function initTopo(condIdx)
        
        
        condData(condIdx).topoH = plotTopo(squeeze( condData(condIdx).Amp(:, condData(condIdx).iFr)),cfg.layout);
        colormap(condData(condIdx).topoAx,hot);
        
        colorbarH = colorbar('peer',condData(condIdx).topoAx,'WestOutside');
        
        colorbarH.Label.String = 'Microvolts';
        
        cpos = colorbarH.Position;
        cpos(3) = 0.5*cpos(3);
        colorbarH.Position = cpos;
        
        set(gcf,'KeyPressFcn',@keyInput)
        
        
        axis off;
        condData(condIdx).elecVerts = get(condData(condIdx).topoH,'Vertices');
        hold on;
        condData(condIdx).markH = plot( condData(condIdx).elecVerts( condData(condIdx).iElec,1), ...
            condData(condIdx).elecVerts( condData(condIdx).iElec,2),...
            'ko','markersize',10,'linewidth',2);
                
        
        set(condData(condIdx).topoH,'ButtonDownFcn',{@clickedTopo,condIdx})
        set(condData(condIdx).topoAx,'ButtonDownFcn',{@clickedTopo,condIdx})
       
        
    end

    function drawInfoPane()
    
        %Using a uitable to present these.  Works well enough. But not the
        %most customizable. 
        
%         infoString = sprintf('TEST\ttest','position')
%         uicontrol('Style','text','units','normalized','parent',infoPanelH,...
%             'string',infoString,'position', [.05 .9 .95 .1]);
        
        rowNames = {...
            'Elec',...
            'Freq',...            
            'Amp',...
            'Phase',...
            'pValue',...
            };
        
        infoPaneData = [...
            condData(1).iElec condData(2).iElec;...
            condData(1).freq(condData(1).iFr) condData(2).freq(condData(2).iFr);...
            condData(1).Amp(condData(1).iElec,condData(1).iFr) condData(2).Amp(condData(2).iElec,condData(2).iFr);... 
            0 0;
            condData(1).tCircPval(condData(1).iElec,condData(1).iFr) condData(2).tCircPval(condData(2).iElec,condData(2).iFr);... 
            ];
        
        %infoPaneData(3,:) = round(infoPaneData(3,:),3,'significant');
        
        %The bank format uses 2 decimals of precision. A bit of a kludge to
        %avoid custom rounding in the above data table values. 
        uitable('Data',infoPaneData,'units','normalized','parent',infoPanelH,...
            'position', [.05 0 .95 .9],'ColumnWidth',{40 40},...
            'columnName',{'A','B'},'rowname',rowNames);
        
        
    
    end
       

    function drawSpec(condIdx)
        
        plotData = condData(condIdx);
        
        plotData.sigFreqs = plotData.tCircPval<=configOptions.pValThresh;
        
        axes(plotData.specAx)
        [plotData.barH plotData.sigH] = pdSpecPlot(plotData.freq,plotData.Amp(plotData.iElec,:),...
            plotData.sigFreqs(plotData.iElec,:));
        title(plotData.specAx,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz'])
        
        %Set up the function to call when the plots are clicked on
        set(plotData.barH,'ButtonDownFcn',{@clickedSpec, condIdx})
        set(plotData.sigH,'ButtonDownFcn',{@clickedSpec, condIdx})
        set(plotData.specAx,'ButtonDownFcn',{@clickedSpec, condIdx})
         drawInfoPane()
    end


    function clickedSpec(hObject,callbackdata,condIdx)
        
        %Get the data to plot
        plotData = condData(condIdx);
         
       tCurPoint = get(plotData.specAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance plotData.iFr] = min(abs(plotData.freq-tCurPoint(1,1)));
        
        if plotData.iFr>=1 && plotData.iFr<=length(plotData.freq)

            set(plotData.topoH,'facevertexCData',plotData.Amp(:,plotData.iFr))
            caxis(plotData.topoAx,[0 max(abs(plotData.Amp(:,plotData.iFr)))])
            title(plotData.specAx,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz'])            
            colormap(plotData.topoAx,hot);            
            drawPhase(plotData);
            condData(condIdx).iFr = plotData.iFr;
             drawInfoPane()
        else
            disp('clicked out of bounds')
        end
    end

    function clickedTopo(hObject,callbackdata,condIdx)
        
        plotData = condData(condIdx);
        
        tCurPoint = get(plotData.topoAx,'CurrentPoint');
        
        dist = bsxfun(@minus,tCurPoint(1,:),plotData.elecVerts);
        
        dist = sqrt(sum(dist.^2,2));
        
        %Get the index to the nearest clicked electrode
        [distance iE] = min(dist);
               
        if iE>=1 && iE<=size(plotData.elecVerts,1),
            
            condData(condIdx).iElec = iE;
            iElec = iE;
            delete(condData(condIdx).markH);
            condData(condIdx).markH = plot(plotData.elecVerts(iElec,1),plotData.elecVerts(iElec,2),'ko','markersize',15,'linewidth',2);
            
            title(condData(condIdx).topoAx,[num2str(iElec) ': ' cfg.layout.label{iElec}]);
            drawSpec(condIdx);
            drawWave(condIdx);
            drawPhase(condIdx);
            drawInfoPane();
        end
        
    end

 function drawWave(condIdx)
     
        %Get the data to plot
        plotData = condData(condIdx);
        
        axes(condData(condIdx).waveAx)
        delete(condData(condIdx).butterflyH);
        condData(condIdx).butterflyH = plot(condData(condIdx).waveAx,condData(condIdx).Time,condData(condIdx).Wave','-','color',[.5 .5 .5],'linewidth',.1);
        hold on;       
        delete(condData(condIdx).selectedLineH);
        condData(condIdx).selectedLineH = plot(condData(condIdx).waveAx,condData(condIdx).Time,...
            condData(condIdx).Wave(condData(condIdx).iElec,:),'color',condData(condIdx).color,'linewidth',2);  
        
        yLims = 1.1*[-max(abs(condData(condIdx).Wave(:)))-eps max(abs(condData(condIdx).Wave(:)))+eps];
        
        axis(condData(condIdx).waveAx,[0 condData(condIdx).Time(end) yLims])
        
        [axLim] = axis(condData(condIdx).waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        delete(condData(condIdx).timeLine);
        condData(condIdx).timeLine = line([condData(condIdx).Time(condData(condIdx).iT) condData(condIdx).Time(condData(condIdx).iT)],[yLo yHi],...
            'linewidth',2,'buttondownFcn',{@clickedWave,condIdx});
        
        ylabel('uV');
        xlabel('Time (ms)');
        %Set up the function to call when the plots are clicked on
        set(condData(condIdx).waveAx,'ButtonDownFcn',{@clickedWave,condIdx})
        set(condData(condIdx).selectedLineH,'ButtonDownFcn',{@clickedWave,condIdx})
        set(condData(condIdx).butterflyH,'ButtonDownFcn',{@clickedWave,condIdx})
        
        
    end

    function clickedWave(hObject,callbackdata,condIdx)
        
        %Get the data to plot
        plotData = condData(condIdx);
        
        tCurPoint = get(plotData.waveAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance iT] = min(abs(plotData.Time-tCurPoint(1,1)));
        
        [axLim] = axis(plotData.waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        
        axes(plotData.waveAx);
        if iT>=1 && iT<=size(plotData.Amp,2)
            
            set(plotData.topoH,'facevertexCData',plotData.Wave(:,iT))
            %           caxis(topoAx,[-max(abs(data(iT,:))) max(abs(data(iT,:)))])
            
            colormap(plotData.topoAx,jmaColors('arizona'));
            caxis(plotData.topoAx,[-max(abs(plotData.Wave(:))) max(abs(plotData.Wave(:)))]);            
            
            set(plotData.timeLine,'XData',[plotData.Time(iT) plotData.Time(iT)],'YData',[yLo yHi]);
            title(plotData.waveAx,['Time: ' num2str(plotData.Time(iT),4) ' ms']);
            
        else
            disp('clicked out of bounds')
        end
    end


    function drawPhase(condIdx)
        %Temp disable while working on other plots
        return;
        axes(phasorAx);
        delete(allchild(phasorAx)); %pdPhasePlot uses weird plotting functions so all it's objects need to be cleared to delete the scale. 
        
        for iCond = 1:2,
        iElec = condData(iCond).iElec;
        iFr = lastFr;
        phaseDataToPlot(iCond) = complex( condData(iCond).Cos(iElec,iFr), condData(iCond).Sin(iElec,iFr));
        end
        
        
        steadyState.tCircStdErr(iElec,iFr)
        
        pdPhasePlot( phaseDataToPlot,steadyState.tCircStdErr(iElec,iFr));
        
    end


    function keyInput(src,evnt)
    
        %disabled while working on other plots. 
        return;
        switch(lower(evnt.Key))
            case 'leftarrow'
                iFr = max(iFr-1,1);
            case 'rightarrow'
                iFr= min(iFr+1,length(freq));
                
                
        end
        
        
          set(topoH,'facevertexCData',pdAmp(:,iFr))
          caxis(topoAx,[0 max(abs(pdAmp(:,iFr)))])
          title(specAx,['Frequency: ' num2str(freq(iFr)) ' Hz'])
            
        
    end


    function selCond(hObject,callbackdata,condIdx)
        
        selectedCond = get(hObject,'Value');
        condData(condIdx) = updateStruct(condData(condIdx),steadyState(selectedCond));
        condData(condIdx).selectedCond = selectedCond;
        
        if condIdx ==1
            configOptions.condIdxA = selectedCond;
        else
            configOptions.condIdxB = selectedCond;
        end
        
        plotData = condData(condIdx);
        set(plotData.topoH,'facevertexCData',plotData.Amp(:,plotData.iFr))
        caxis(plotData.topoAx,[0 max(abs(plotData.Amp(:,plotData.iFr)))])
        title(plotData.specAx,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz'])
        colormap(plotData.topoAx,hot);
        
        
        
        %         cla(condData(condIdx).topoAx);
        %         initTopo(condIdx);
        
        cla(condData(condIdx).specAx);
        drawSpec(condIdx);
        
        cla(condData(condIdx).waveAx);
        drawWave(condIdx);
        
        drawPhase(condIdx);
        
    end

    function selPVal(hObject,callbackdata)
        
        %Grab the selected menu item.
        pValIdx = get(hObject,'Value');
        popupList = get(hObject,'String');
        selectedPVal = popupList{pValIdx};
        
        if strcmpi(selectedPVal,'none')
            configOptions.pValThresh = Inf;
        else
            configOptions.pValThresh = str2double(selectedPVal);
        end
        
        
        drawSpec(1);
        drawSpec(2);
    end

    function selFilter(hObject,callbackdata,condIdx)
        
        %Grab the selected menu item.
        filtIdx = get(hObject,'Value');
        filtList = get(hObject,'String');
        filterName = filtList{filtIdx};
        
        filterWave(filterName,condIdx);
        
        drawSpec(condIdx);        
        drawWave(1);
        drawWave(2);
        
    end

    function filterWave(filterName,condIdx)

        
        %Determine what filter values to keep. 
        if strcmpi(filterName,'none')
            filtIdx = 2:condData(condIdx).nFr; %Do not include DC for plot. 
        else
            filtIdx = determineFilterIndices(filterName,condData(condIdx).freq,condData(condIdx).i1F1);
        end
        filtIdx
        cfg.keepFreqIdx = filtIdx;
        condData(condIdx).Wave = filterSteadyState(cfg,steadyState(condData(condIdx).selectedCond));
        

        
    end



end
