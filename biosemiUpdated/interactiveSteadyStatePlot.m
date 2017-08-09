function [] = interactiveSteadyStatePlot(cfg, steadyState)
%function [] = interactiveTopoPlot(cfg,steadyState)
%
% helpful help here
%


%Setup default condition choices.
configOptions.selCondIdx(1) = 1;
configOptions.selCondIdx(2) = 1;
if length(steadyState)>=2, %If 2 or more conditions default to showing cond 2. 
    configOptions.selCondIdx(2) = 2;
end

configOptions.pValThresh = Inf;

%Options for different automatic yscaling methods:

configOptions.yScale = 'global'; %Lock y-scales by the largest in any condition
%configOptions.yScale = 'local'; %Y-scales maximize each independent. 
%configOptions.yScale = 'byCondition'; %Lock y-scales for each condition indepenedtly

configOptions.underlay = 'none'; %What data should be shown undet the waveform
configOptions.compareWave = false; %Should we overlay wave forms?

%TODO: Clarify and check this code. 
% prepare the layout, this should be done only once
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg);



%
%configOptions - A structure containing the different options to use for
%                displaying data
%condData
%


%Using condData to hold th 2 conditions to plot. 

condData(1) = initPlotData(steadyState(configOptions.selCondIdx(1)));
condData(2) = initPlotData(steadyState(configOptions.selCondIdx(2)));
condData(1).selectedCond = configOptions.selCondIdx(1);
condData(2).selectedCond = configOptions.selCondIdx(2);


%Sets the colors for the two conditions. 
condData(1).color = [1 0 0];
condData(2).color = [0 0 1];


iElec = 1;
setPlotGlobalLimits();

%Initialize figure;
figH= figure('units','normalized','outerposition',[0 1 .8 1]); %Render a new figure.

%Set the color order for this figure to mimic the condition colors.
set(figH,'DefaultAxesColorOrder',[condData(1).color; condData(2).color]);



%%%%% Setup plot option GUI controls
% Create pop-up menu

%Quickie Just label conditions with numbers:
conditionList = cellstr(num2str([1:length(steadyState)]'));

pValList = {'None','0.05','0.01','1e-3','1e-6'};
filterList = {'None','nF1','nF1Odd','nF1Even','nF1Low49'};
    
uicontrol('Style', 'text','units','normalized',...
    'String','Red',...
    'Position', [.9 .975 .05 .02])

uicontrol('Style', 'text','units','normalized',...
    'String','Blue',...
    'Position', [.95 .975 .05 .02])

uicontrol('Style', 'text','units','normalized',...
    'String','Condition: ',...
    'Position', [.85 .94 .05 .02])

selCondAPopupH = uicontrol('Style', 'popup','units','normalized',...
    'String',conditionList,'Value',configOptions.selCondIdx(1),...
    'Position', [.9 .95 .05 .01],...
    'Callback', {@selCond,1});

selCondBPopupH = uicontrol('Style', 'popup','units','normalized',...
    'String',conditionList,'Value',configOptions.selCondIdx(2),...
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

uicontrol('Style', 'text','units','normalized',...
    'String','Wave Scale:',...
    'Position', [.81 .8 .06 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',{'Global Max','Local Max','Per Condition Global'},'Value',1,...
    'Position', [.87 .8 .1 .02],'fontsize',8,...
    'Callback', @selYlim);

uicontrol('Style', 'text','units','normalized',...
    'String','Underlay:',...
    'Position', [.81 .75 .06 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',{'None','Butterfly'},'Value',1,...
    'Position', [.87 .75 .1 .02],'fontsize',8,...
    'Callback', @selUnderlay);

uicontrol('Style', 'text','units','normalized',...
    'String','Compare:',...
    'Position', [.81 .71 .06 .02]);

uicontrol('Style', 'checkbox','units','normalized',...
    'String','','Value',false,...
    'Position', [.87 .71 .1 .02],'fontsize',8,...
    'Callback', @selCompare);


infoPanelH = uipanel('Title','Info','FontSize',12,...
             'BackgroundColor','white',...
             'Position',[.82 .4 .16 .4],'visible','off');


%%%%%%% Setup data plot axes.

%Creating axes here to setup layout for where things are plotted. 
%Default Topography is spec. 

plotXStart = .35;

%topoAx = subplot(10,1,1:8);
condData(1).topoAx = axes('Parent',figH,'Position',[0.05 .6 .3 .3]);
initTopo(1);

condData(2).topoAx = axes('Parent',figH,'Position',[0.05 .25 .3 .3]);
initTopo(2);

%Setup the frequency domain plot
%For condition A
condData(1).specAx = axes('Parent',figH,'Position',[plotXStart .82 .4 .12]);
drawSpec(1);

%For condition B
condData(2).specAx = axes('Parent',figH,'Position',[plotXStart .57 .4 .12]);
drawSpec(2);

%Setup the time domain plot
condData(1).waveAx = axes('Parent',figH,'Position',[plotXStart .29 .4 .15]);
drawWave(1);

condData(2).waveAx = axes('Parent',figH,'Position',[plotXStart .08 .4 .15]);
drawWave(2);



%Setup complex phasor plot
phasorAx = axes('Parent',figH,'Position',[0.73 .1 .3 .3]);
drawPhase();


%%%% Setup info pane information
drawInfoPane()


%set(topoAx,'ButtonDownFcn',@specUpdate)


    function plotData = initPlotData(steadyState)
    %This function does the setup needed to take the steadyState fields and
    %make them into values needed for plotting. It also adds fields that 
    %are useful for plotting data. 
        
    plotData = steadyState;
    
    %!!!!!! This is a really, really stupid line.  Just relying on stupidy
    %being obvious on the plot if units are off by x1000
    plotData.time = 1000*steadyState.time; %TODO: Make time units more explicit.        
    
    if ~isfield(steadyState,'pval')
        plotData.sigFreqs = [];
    else
        plotData.sigFreqs = steadyState.pval<=configOptions.pValThresh;
    end
    
    plotData.filterName = 'None';
    plotData.activeFreq = true(size(plotData.amp)); 
    plotData.activeFreq(:,1) = false;%Turn OFF DC by default;
    
    %Set default plot options        
    plotData.iElec = 1;
    plotData.iFr = 1;
    plotData.iT = 1;
    
    
    
    %Handles for plot elements
    plotData.timeLine = [];
    plotData.butterflyH = [];
    plotData.selectedLineH = [];
    plotData.overlayLineH = [];
    
    end


    function setPlotGlobalLimits()
        %This function sets plot limits when global scaling is used
        
        for iCond = 1:2,
            condIdx = configOptions.selCondIdx(iCond);
            waveMax(iCond) = max(abs(steadyState(condIdx).wave(:)));
            specMax(iCond) = max(abs(steadyState(condIdx).amp(:)));
        end
        
        %If we want global max, replace values with global max
        if strcmpi(configOptions.yScale,'global')            
            waveMax(1:2) = max(waveMax);
            specMax(1:2) = max(specMax);
        end
        
        for iCond = 1:2,
        %Now loop again and set the values.
            %Have to add/subtract eps for cases in which plot data is all flat
            %(e.g. ref channel
        condData(iCond).globalWaveYLims  = 1.1*[-waveMax(iCond)-eps waveMax(iCond)+eps];
        condData(iCond).globalSpecYLims  = 1.1*[0 specMax(iCond)+eps];
        end
        
    end

    function initTopo(condIdx)
        
        
        condData(condIdx).topoH = plotTopo(squeeze( condData(condIdx).amp(:, condData(condIdx).iFr)),cfg.layout);
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
            condData(1).amp(condData(1).iElec,condData(1).iFr) condData(2).amp(condData(2).iElec,condData(2).iFr);... 
            0 0;
            condData(1).pval(condData(1).iElec,condData(1).iFr) condData(2).pval(condData(2).iElec,condData(2).iFr);... 
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
        
   
        
        sigFreqs =  plotData.pval(plotData.iElec,:)<=configOptions.pValThresh;
        axes(plotData.specAx)
        [plotData.barH plotData.sigH] = pdSpecPlot(plotData.freq,plotData.amp(plotData.iElec,:),...
            sigFreqs);
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

            
            condData(condIdx).topoMode = 'amp';
            title(plotData.specAx,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz']);            

            condData(condIdx).iFr = plotData.iFr;
                        
            drawPhase();            
            drawTopo(condIdx);
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
            refreshPlots();
        end
        
    end

    function drawTopo(condIdx);
        %Draw top
        plotData = condData(condIdx);
        
        
        switch lower(condData(condIdx).topoMode)
            case 'amp'
        
                maxVal = max(abs(plotData.amp(:,plotData.iFr)));
                %                 if strcmpi(configOptions.yScale,'local') %If local use the locally deterrmined maxval
                %                     yLims = 1.1*[0 maxVal+eps]; %Set the yLimits
                %                 else
                %                     yLims = plotData.globalSpecYLims;
                %                 end
                    
                %Always using local scaling for AMP topo plots, makes most
                %sense in most use cases.  So disabling user selectable
                %scaling for now. 
                yLims = 1.1*[0 maxVal+eps]; %Set the yLimits

                set(plotData.topoH,'facevertexCData',plotData.amp(:,plotData.iFr));
                caxis(plotData.topoAx,yLims);
                colormap(plotData.topoAx,hot);
            
            case 'wave'
                
                maxVal = max(abs(plotData.wave(:)));
                if strcmpi(configOptions.yScale,'local') %If local use the locally deterrmined maxval
                    yLims = 1.1*[-maxVal-eps maxVal+eps]; %Set the yLimits
                else
                    yLims = plotData.globalWaveYLims;
                end
        
                set(plotData.topoH,'facevertexCData',plotData.wave(:,plotData.iT))                                
                caxis(plotData.topoAx,yLims);
                colormap(plotData.topoAx,jmaColors('arizona'));
        end
        
    end

    function drawWave(condIdx)
     
        %Get the data to plot
        plotData = condData(condIdx);
        condData(condIdx).iT = min(condData(condIdx).nt,condData(condIdx).iT);
        
        axes(condData(condIdx).waveAx)
        delete(condData(condIdx).butterflyH);
        delete(condData(condIdx).selectedLineH);
        delete(condData(condIdx).overlayLineH);
       
        selectedWave = condData(condIdx).wave(condData(condIdx).iElec,:);
        maxVal = max(selectedWave);
        
        if strcmp(configOptions.underlay,'butterfly')
            condData(condIdx).butterflyH = plot(condData(condIdx).waveAx,condData(condIdx).time,condData(condIdx).wave','-','color',[.5 .5 .5],'linewidth',.1);
            set(condData(condIdx).butterflyH,'ButtonDownFcn',{@clickedWave,condIdx});
            maxVal = max(abs(condData(condIdx).wave(:)));
        end
        
        
        condData(condIdx).selectedLineH = plot(condData(condIdx).waveAx,condData(condIdx).time,...
           selectedWave,'color',condData(condIdx).color,'linewidth',2);  
        hold on;
        
        if configOptions.compareWave
            overlayIdx = -condIdx+3; %Tricky way to change between 2 and 1 
            overlayWave = condData(overlayIdx).wave(condData(overlayIdx).iElec,:);
            condData(condIdx).overlayLineH =plot(condData(condIdx).waveAx,condData(overlayIdx).time,...
                overlayWave,'color',condData(overlayIdx).color,'linewidth',2);
        end
        
        if strcmpi(configOptions.yScale,'local') %If local use the locally deterrmined maxval
            yLims = 1.1*[-maxVal-eps maxVal+eps]; %Set the yLimits
        else
            yLims = plotData.globalWaveYLims;
        end
        
        axis(condData(condIdx).waveAx,[0 condData(condIdx).time(end) yLims])
        
        [axLim] = axis(condData(condIdx).waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        delete(condData(condIdx).timeLine);
        condData(condIdx).timeLine = line([condData(condIdx).time(condData(condIdx).iT) condData(condIdx).time(condData(condIdx).iT)],[yLo yHi],...
            'linewidth',2,'buttondownFcn',{@clickedWave,condIdx});
        
        ylabel('uV');
        
        %Only plot the xlabel if rendering the bottom plot.  That enables
        %us to make the plots bigger
        xlabel('');
        if condIdx ==2;
        xlabel('Time (ms)');
        end
        %Set up the function to call when the plots are clicked on
        set(condData(condIdx).waveAx,'ButtonDownFcn',{@clickedWave,condIdx})
        set(condData(condIdx).selectedLineH,'ButtonDownFcn',{@clickedWave,condIdx})
        
        
    end

    function clickedWave(hObject,callbackdata,condIdx)
        
        %Get the data to plot
        plotData = condData(condIdx);
        
        tCurPoint = get(plotData.waveAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance iT] = min(abs(plotData.time-tCurPoint(1,1)));
        
        [axLim] = axis(plotData.waveAx);
        yLo = axLim(3);
        yHi = axLim(4);
        
        axes(plotData.waveAx);
        if iT>=1 && iT<=size(plotData.wave,2)
            
            
            set(plotData.timeLine,'XData',[plotData.time(iT) plotData.time(iT)],'YData',[yLo yHi]);
            title(plotData.waveAx,['Time: ' num2str(plotData.time(iT),4) ' ms']);
            
            condData(condIdx).topoMode = 'wave';
            condData(condIdx).iT = iT;
            drawTopo(condIdx);
                        
        else
            disp('clicked out of bounds')
        end
    end


    function drawPhase(condIdx)
        
        
        axes(phasorAx);
        delete(allchild(phasorAx)); %pdPhasePlot uses weird plotting functions so all it's objects need to be cleared to delete the scale. 
        
        %First lets make the data complex.
        for iCond = 1:2,
            iElec = condData(iCond).iElec;
            iFr   = condData(iCond).iFr;
            phaseDataToPlot(iCond) = complex( condData(iCond).cos(iElec,iFr), condData(iCond).sin(iElec,iFr));
            
            %TODO: Fix this! this is an underestimate of the 95% confidence
            %intervals!
            confRadius(iCond) = condData(iCond).confradius(iElec,iFr);
        end
                      
        
        pdPhasePlot( phaseDataToPlot,confRadius);
        
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

    %Callback when user selects condition from dropdown. 
    function selCond(hObject,callbackdata,condIdx)
        
        %Boilerplate code finding what was selected.
        selectedCond = get(hObject,'Value');
        
        %Update struct will take the data from steadyState and update the
        %condData. 
        condData(condIdx) = updateStruct(condData(condIdx),steadyState(selectedCond));
        condData(condIdx).selectedCond = selectedCond;
        
        %Set the selected condtion. 
        configOptions.selCondIdx(condIdx) = selectedCond;
        
        plotData = condData(condIdx);
        set(plotData.topoH,'facevertexCData',plotData.amp(:,plotData.iFr))
        caxis(plotData.topoAx,[0 max(abs(plotData.amp(:,plotData.iFr)))])
        title(plotData.specAx,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz'])
        colormap(plotData.topoAx,hot);
        
        filterWave(condIdx)
        refreshPlots();
        
        %         cla(condData(condIdx).topoAx);
        %         initTopo(condIdx);
        
%         cla(condData(condIdx).specAx);
%         drawSpec(condIdx);
%         
%         cla(condData(condIdx).waveAx);
%         drawWave(condIdx);
%         
%         drawPhase(condIdx);
        
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
        
        filterWave(1);
        filterWave(2);
        
        drawSpec(1);
        drawSpec(2);
        drawWave(1);
        drawWave(2);
    end

    function selFilter(hObject,callbackdata,condIdx)
        
        %Grab the selected menu item.
        filtIdx = get(hObject,'Value');
        filtList = get(hObject,'String');
        filterName = filtList{filtIdx};
        
        condData(condIdx).filterName = filterName;
        filterWave(condIdx);
        
        drawSpec(condIdx);        
        drawWave(1);
        drawWave(2);
        
    end


    function selYlim(hObject,callbackdata)
        %Grab the selected menu item.
        idx = get(hObject,'Value');
        list = get(hObject,'String');
        limName = list{idx};
        
        switch limName
            
            case 'Global Max',
                configOptions.yScale = 'global'; %Lock y-scales by the largest in any condition
            case 'Local Max',
                configOptions.yScale = 'local'; %Y-scales maximize each independent.
            case 'Per Condition Global',
                configOptions.yScale = 'byCondition'; %Lock y-scales for each condition indepenedtly
        end
        
        setPlotGlobalLimits();
        refreshPlots();
    end


    function selUnderlay(hObject,callbackdata)
        %Select what underlay data to show. 
        
        %Grab the selected menu item.
        idx = get(hObject,'Value');
        list = get(hObject,'String');
        choice = list{idx};
        
            
        configOptions.underlay = lower(choice);
        
        drawWave(1);
        drawWave(2);
    end


    function selCompare(hObject,callbackdata,condIdx)
        
        configOptions.compareWave = get(hObject,'Value');
        drawWave(1);
        drawWave(2);
        
    end

    function refreshPlots()
        %Refresh all the data plots.
        
        drawSpec(1);
        drawSpec(2);
        drawWave(1);
        drawWave(2);
        drawPhase();
    end


    function filterWave(condIdx)
        
        
        %Create a logical matrix selecting significant components.
        sigFreqs =  condData(condIdx).pval<=configOptions.pValThresh;
        
        
        %Determine what filter values to keep.
        %This does different things if 'none' is chosen.
        if strcmpi(condData(condIdx).filterName,'none')
            filtIdx = 2:condData(condIdx).nfr; %Do not include DC.
            condData(condIdx).wave = steadyState(condData(condIdx).selectedCond).wave;
        else
            filtIdx = determineFilterIndices(condData(condIdx).filterName,...
                condData(condIdx).freq, condData(condIdx).i1f1);
        end
        
        %Create a logical matrix selecting frequency components.
        filtMat = false(size(condData(condIdx).amp));
        filtMat(:,filtIdx) = true;
        
        
        %Combine the filter and sig vaules with a logical AND.
        condData(condIdx).activeFreq = (filtMat.*sigFreqs)>=1; %Store the filtered coefficients for the spec plot
        

        cfg.activeFreq =  condData(condIdx).activeFreq;
 
        %if 'none' is chosen don't apply any filter
        if strcmpi(condData(condIdx).filterName,'none')
            
            condData(condIdx).wave = steadyState(condData(condIdx).selectedCond).wave;
        else
            condData(condIdx).wave = filterSteadyState(cfg,steadyState(condData(condIdx).selectedCond));
        end
        
        

        
    end



end
