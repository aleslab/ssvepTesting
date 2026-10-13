function [] = interactiveSteadyStatePlot2(cfg, steadyState)
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
configOptions.lowpassHz = 49;
configOptions.applyLowpass = false;
configOptions.showPTopo = false;

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


setPlotGlobalLimits();

%Initialize figure;
figH= figure('WindowStyle','normal','units','normalized','Position',[.02 .05 .94 .8]); %Render a new figure.

% Let MATLAB finish sizing the client area before placing normalized controls.
drawnow;

%Set the color order for this figure to mimic the condition colors.
set(figH,'DefaultAxesColorOrder',[condData(1).color; condData(2).color]);



%%%%% Setup plot option GUI controls
% Create pop-up menu
if isfield(steadyState,'condLabel')
    conditionList = {steadyState.condLabel}';
else
    %Quickie Just label conditions with numbers:
    conditionList = cellstr(num2str([1:length(steadyState)]'));
end

pValList = {'None','0.05','0.01','1e-3','1e-6'};
filterList = {'None','nF1','nF1Odd','nF1Even','nF1Clean', ...
    'nF2','nF2Odd','nF2Even','nF2Clean'};
    
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
    'Position', [.9 .93 .05 .04],...
    'Callback', {@selCond,1});

selCondBPopupH = uicontrol('Style', 'popup','units','normalized',...
    'String',conditionList,'Value',configOptions.selCondIdx(2),...
    'Position', [.95 .93 .05 .04],...
    'Callback', {@selCond,2});


uicontrol('Style', 'text','units','normalized',...
    'String','p Threshold:',...
    'Position', [.82 .9 .08 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',pValList,'Value',1,...
    'Position', [.9 .9 .1 .04],...
    'Callback', @selPVal);

uicontrol('Style', 'text','units','normalized',...
    'String','Time Filter:',...
    'Position', [.81 .85 .06 .02]);

filterPopupH(1) = uicontrol('Style', 'popup','units','normalized',...
    'String',filterList,'Value',1,...
    'Position', [.87 .85 .07 .04],'fontsize',8,...
    'Callback', {@selFilter,1});

filterPopupH(2) = uicontrol('Style', 'popup','units','normalized',...
    'String',filterList,'Value',1,...
    'Position', [.93 .85 .07 .04],'fontsize',8,...
    'Callback', {@selFilter,2});

uicontrol('Style', 'text','units','normalized',...
    'String','Wave Scale:',...
    'Position', [.81 .8 .06 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',{'Global Max','Local Max','Per Condition Global'},'Value',1,...
    'Position', [.87 .8 .1 .04],'fontsize',8,...
    'Callback', @selYlim);

uicontrol('Style', 'text','units','normalized',...
    'String','Underlay:',...
    'Position', [.81 .75 .06 .02]);

uicontrol('Style', 'popup','units','normalized',...
    'String',{'None','Butterfly'},'Value',1,...
    'Position', [.87 .75 .1 .04],'fontsize',8,...
    'Callback', @selUnderlay);

uicontrol('Style', 'text','units','normalized',...
    'String','Compare:',...
    'Position', [.81 .71 .06 .02]);

uicontrol('Style', 'checkbox','units','normalized',...
    'String','','Value',false,...
    'Position', [.87 .71 .1 .02],'fontsize',8,...
    'Callback', @selCompare);


windowLabels = {'Red window:','Blue window:'};
for k = 1:2
    uicontrol('Style','text','units','normalized','String',windowLabels{k}, ...
        'Position',[.81 .65-(k-1)*.05 .08 .025]);
    wavePopupH(k) = uicontrol('Style','popup','units','normalized', ...
        'String',condData(k).wavelabels,'Value',1, ...
        'Position',[.89 .65-(k-1)*.05 .11 .04], ...
        'Callback',{@selWave,k});
end

uicontrol('Style','text','units','normalized','String','Low-pass (Hz):', ...
    'Position',[.81 .54 .09 .025]);
lowpassEditH = uicontrol('Style','edit','units','normalized', ...
    'String',num2str(configOptions.lowpassHz),'Position',[.90 .54 .045 .035], ...
    'FontSize',9,'Callback',@setLowpassCutoff);
uicontrol('Style','togglebutton','units','normalized', ...
    'String','Apply','Value',configOptions.applyLowpass, ...
    'Position',[.95 .54 .05 .035],'FontSize',9,'Callback',@toggleLowpass);

uicontrol('Style','togglebutton','units','normalized', ...
    'String','P-value topography','Value',configOptions.showPTopo, ...
    'Position',[.84 .48 .15 .035],'FontSize',9,'Callback',@togglePTopo, ...
    'TooltipString','Show p-values at the selected frequency; grey means p exceeds the threshold.');

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
    
    % New analyses explicitly provide ms. Older analysis output has dtms;
    % legacy imported data without this metadata retains the seconds convention.
    if ~isfield(plotData,'waveforms') || isempty(plotData.waveforms)
        if isfield(plotData,'dtms') || ...
                (isfield(plotData,'timeunit') && strcmp(plotData.timeunit,'ms'))
            timeMs = plotData.time;
        else
            timeMs = 1000*plotData.time;
        end
        plotData.waveforms = {plotData.wave};
        plotData.wavetimes = {timeMs};
        plotData.wavelabels = {'Default'};
    end
    plotData.waveIdx = 1;
    plotData.wave = plotData.waveforms{1};
    plotData.time = plotData.wavetimes{1};
    plotData.nt = size(plotData.wave,2);
    plotData.topoMode = 'amp';

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
        %Waveform limits only; spectra always scale to their displayed data.
        
        for iCond = 1:2,
            waveMax(iCond) = max(abs(condData(iCond).wave(:)));
        end
        
        %If we want global max, replace values with global max
        if strcmpi(configOptions.yScale,'global')            
            waveMax(1:2) = max(waveMax);
        end
        
        for iCond = 1:2,
        %Now loop again and set the values.
            %Have to add/subtract eps for cases in which plot data is all flat
            %(e.g. ref channel
        condData(iCond).globalWaveYLims  = 1.1*[-waveMax(iCond)-eps waveMax(iCond)+eps];
        end
        
    end

    function initTopo(condIdx)
        
        
        condData(condIdx).topoH = plotTopo(squeeze( condData(condIdx).amp(:, condData(condIdx).iFr)),cfg.layout);
        colormap(condData(condIdx).topoAx,hot);
        
        colorbarH = colorbar('peer',condData(condIdx).topoAx,'WestOutside');
        
        colorbarH.Label.String = 'Microvolts';
        condData(condIdx).topoColorbarH = colorbarH;
        
        cpos = colorbarH.Position;
        cpos(3) = 0.5*cpos(3);
        colorbarH.Position = cpos;
        
        set(gcf,'KeyPressFcn',@keyInput)
        
        
        axis off;
        condData(condIdx).elecVerts = get(condData(condIdx).topoH,'Vertices');
        hold(condData(condIdx).topoAx,'on');
        condData(condIdx).markH = plot(condData(condIdx).topoAx, condData(condIdx).elecVerts( condData(condIdx).iElec,1), ...
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
        ax = plotData.specAx;
        % Explicitly replace all spectrum objects, including old legend keys.
        % Never rely on the current axes or inherited hold state.
        legend(ax,'off');
        cla(ax);
        hold(ax,'on');
        amps = plotData.amp(plotData.iElec,:);
        sigFreqs = significantBins(plotData);
        selected = plotData.activeFreq(plotData.iElec,:) & sigFreqs(plotData.iElec,:);
        selected(1) = false;
        displayAmps = amps;
        displayAmps(1) = 0; % DC remains available for selection, not y scaling.
        % One bar series prevents grouped/overlaid bars from shifting or
        % obscuring bins when MATLAB redraws the spectrum.
        condData(condIdx).barH = bar(ax,plotData.freq,displayAmps,.8, ...
            'FaceColor','flat','EdgeColor','none');
        colourSpectrum(condIdx,selected);
        upperHz = min(49,max(plotData.freq));
        xlim(ax,[0 upperHz]);
        % Spectrum scaling is independent of the Wave Scale control.
        shown = plotData.freq > 0 & plotData.freq <= upperHz;
        ylim(ax,[0 1.1*(max(amps(shown))+eps)]);
        title(ax,['Frequency: ' num2str(plotData.freq(plotData.iFr)) ' Hz']);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'uV');
        set(condData(condIdx).barH,'ButtonDownFcn',{@clickedSpec,condIdx});
        set(ax,'ButtonDownFcn',{@clickedSpec,condIdx});
        hold(ax,'off');
    end

    function mask = significantBins(plotData)
        mask = true(size(plotData.amp));
        if isfinite(configOptions.pValThresh)
            if ~isfield(plotData,'pval') || isempty(plotData.pval)
                error('interactiveSteadyStatePlot2:MissingPValues', ...
                    'A p-value threshold requires pval in the analysis output.');
            end
            mask = plotData.pval <= configOptions.pValThresh;
        end
    end

    function clickedSpec(hObject,callbackdata,condIdx)
        plotData = condData(condIdx);
        point = get(plotData.specAx,'CurrentPoint');
        [~,iFr] = min(abs(plotData.freq-point(1,1)));
        setFrequencyIndex(iFr);
        for side = 1:2
            condData(side).topoMode = 'amp';
            title(condData(side).specAx, ...
                ['Frequency: ' num2str(condData(side).freq(condData(side).iFr)) ' Hz']);
            drawTopo(side);
        end
        drawPhase();
        drawInfoPane();
    end

    function setFrequencyIndex(iFr)
        % Conditions share the experiment's frequency grid. Keep one bin
        % selected for both spectra, topographies and the phasor comparison.
        iFr = min([iFr numel(condData(1).freq) numel(condData(2).freq)]);
        condData(1).iFr = iFr;
        condData(2).iFr = iFr;
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
            delete(condData(condIdx).markH);
            condData(condIdx).markH = plot(plotData.topoAx,plotData.elecVerts(iE,1),plotData.elecVerts(iE,2),'ko','markersize',15,'linewidth',2);
            
            title(condData(condIdx).topoAx,[num2str(iE) ': ' cfg.layout.label{iE}]);
            refreshPlots();
        end
        
    end

    function drawTopo(condIdx);
        %Draw top
        plotData = condData(condIdx);
        
        
        if configOptions.showPTopo
            drawPTopo(condIdx);
            return
        end
        plotData.topoColorbarH.Label.String = 'Microvolts';
        title(plotData.topoAx,sprintf('%d: %s',plotData.iElec, ...
            cfg.layout.label{plotData.iElec}));
        switch lower(condData(condIdx).topoMode)
            case 'amp'
        
                maxVal = max(abs(plotData.amp(:,plotData.iFr)));
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

    function togglePTopo(hObject,~)
        configOptions.showPTopo = logical(get(hObject,'Value'));
        drawTopo(1);
        drawTopo(2);
    end

    function drawPTopo(condIdx)
        % Use the same significance scale for both selected conditions.
        maxScore = 1;
        for side = 1:2
            p = condData(side).pval(:,condData(side).iFr);
            valid = isfinite(p) & p >= 0 & p <= 1;
            scores = -log10(max(p(valid),realmin));
            maxScore = max([maxScore; scores(:)]);
        end
        plotData = condData(condIdx);
        p = plotData.pval(:,plotData.iFr);
        coloured = isfinite(p) & p >= 0 & p <= 1 & p <= configOptions.pValThresh;
        palette = viridisMap();
        rgb = repmat([.65 .65 .65],numel(p),1);
        scores = -log10(max(p(coloured),realmin));
        colourIdx = 1+round(255*scores/maxScore);
        rgb(coloured,:) = palette(colourIdx,:);
        % Truecolour keeps excluded electrodes grey, outside the colour scale.
        set(plotData.topoH,'FaceVertexCData',rgb);
        colormap(plotData.topoAx,palette);
        caxis(plotData.topoAx,[0 maxScore]);
        plotData.topoColorbarH.Label.String = '-log10(p)';
        title(plotData.topoAx,sprintf('p-values: %.2f Hz',plotData.freq(plotData.iFr)));
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
        maxVal = max(abs(selectedWave));
        
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


    function drawPhase()
        
        
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
        title(phasorAx,'');
        text(phasorAx,0,1.02,sprintf('Phasors: %.2f Hz', ...
            condData(1).freq(condData(1).iFr)), ...
            'Units','normalized','HorizontalAlignment','left', ...
            'VerticalAlignment','bottom','FontSize',8,'FontWeight','normal', ...
            'Clipping','off','HitTest','off');
        
    end


%     function keyInput(src,evnt)
%     
%         %disabled while working on other plots. 
%         return;
%         switch(lower(evnt.Key))
%             case 'leftarrow'
%                 iFr = max(iFr-1,1);
%             case 'rightarrow'
%                 iFr= min(iFr+1,length(freq));
%                 
%                 
%         end
%         
%         
%           set(topoH,'facevertexCData',pdAmp(:,iFr))
%           caxis(topoAx,[0 max(abs(pdAmp(:,iFr)))])
%           title(specAx,['Frequency: ' num2str(freq(iFr)) ' Hz'])
%             
%         
%     end

    %Callback when user selects condition from dropdown. 
    function selCond(hObject,callbackdata,condIdx)
        
        %Boilerplate code finding what was selected.
        selectedCond = get(hObject,'Value');
        
        %Update struct will take the data from steadyState and update the
        %condData. 
        oldLabel = condData(condIdx).wavelabels{condData(condIdx).waveIdx};
        newData = initPlotData(steadyState(selectedCond));
        % Keep graphics handles and filter selection while replacing data.
        fields = fieldnames(steadyState(selectedCond));
        for k = 1:numel(fields)
            condData(condIdx).(fields{k}) = newData.(fields{k});
        end
        condData(condIdx).waveforms = newData.waveforms;
        condData(condIdx).wavetimes = newData.wavetimes;
        condData(condIdx).wavelabels = newData.wavelabels;
        idx = find(strcmp(oldLabel,newData.wavelabels),1);
        if isempty(idx), idx = 1; end
        condData(condIdx).waveIdx = idx;
        setFrequencyIndex(condData(condIdx).iFr);
        condData(condIdx).iElec = min(condData(condIdx).iElec,size(newData.amp,1));
        set(wavePopupH(condIdx),'String',newData.wavelabels,'Value',idx);
        condData(condIdx).selectedCond = selectedCond;
        if strncmpi(condData(condIdx).filterName,'nF2',3) && ...
                (~isfield(newData,'i1f2') || isempty(newData.i1f2))
            condData(condIdx).filterName = 'None';
            set(filterPopupH(condIdx),'Value',1);
            warning('interactiveSteadyStatePlot2:FilterReset', ...
                'This condition has no f2 tag; its time filter has been reset to None.');
        end
        
        %Set the selected condtion. 
        configOptions.selCondIdx(condIdx) = selectedCond;
        
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
        
        refreshPlots();
    end

    function selFilter(hObject,callbackdata,condIdx)
        
        %Grab the selected menu item.
        filtIdx = get(hObject,'Value');
        filtList = get(hObject,'String');
        filterName = filtList{filtIdx};
        
        condData(condIdx).filterName = filterName;
        filterWave(condIdx);
        
        refreshPlots();
        
    end


    function setLowpassCutoff(hObject,~)
        cutoff = str2double(get(hObject,'String'));
        if ~isfinite(cutoff) || cutoff <= 0
            set(lowpassEditH,'String',num2str(configOptions.lowpassHz));
            errordlg('Enter a positive cutoff frequency in Hz.','Low-pass cutoff');
            return
        end
        configOptions.lowpassHz = cutoff;
        filterWave(1);
        filterWave(2);
        refreshPlots();
    end

    function toggleLowpass(hObject,~)
        configOptions.applyLowpass = logical(get(hObject,'Value'));
        filterWave(1);
        filterWave(2);
        refreshPlots();
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
        setPlotGlobalLimits();
        drawTopo(1);
        drawTopo(2);
        
        drawSpec(1);
        drawSpec(2);
        drawWave(1);
        drawWave(2);
        drawPhase();
        drawInfoPane();
    end


    function selWave(hObject,~,condIdx)
        condData(condIdx).waveIdx = get(hObject,'Value');
        filterWave(condIdx);
        refreshPlots();
    end

    function colourSpectrum(condIdx,selected)
        plotData = condData(condIdx);
        % Use exactly the same integer bins as the time filters. Nominal Hz
        % can differ slightly from measured/resampled frequencies.
        n = numel(plotData.freq);
        bins = 0:n-1;
        h1 = false(1,n); h2 = h1; im = h1;
        h1(determineFilterIndices('nF1',plotData.freq,plotData.i1f1)) = true;
        palette = [0 .45 .74; .85 .33 .1; .49 .18 .56];
        labels = {'f1 harmonics'};
        if isfield(plotData,'i1f2') && ~isempty(plotData.i1f2)
            h2(determineFilterIndices('nF2',plotData.freq,plotData.i1f1,plotData.i1f2)) = true;
            commonBin = gcd(plotData.i1f1-1,plotData.i1f2-1);
            im = bins > 0 & mod(bins,commonBin)==0 & ~h1 & ~h2;
            labels = {'f1 (including shared)','f2 only','Intermodulation only'};
        end
        colours = repmat([.6 .6 .6],n,1);
        colours(h1,:) = repmat(palette(1,:),sum(h1),1);
        onlyF2 = h2 & ~h1;
        colours(onlyF2,:) = repmat(palette(2,:),sum(onlyF2),1);
        colours(im,:) = repmat(palette(3,:),sum(im),1);
        colours(~selected,:) = repmat([.85 .85 .85],sum(~selected),1);
        set(plotData.barH,'CData',colours);
        keys = gobjects(1,numel(labels)+1);
        for k = 1:numel(labels)
            keys(k) = plot(plotData.specAx,nan,nan,'s','Color',palette(k,:), ...
                'MarkerFaceColor',palette(k,:),'HitTest','off');
        end
        keys(end) = plot(plotData.specAx,nan,nan,'s','Color',[.6 .6 .6], ...
            'MarkerFaceColor',[.85 .85 .85],'HitTest','off');
        labels{end+1} = 'Outside filter / threshold';
        legend(plotData.specAx,keys,labels,'Location','northeast','FontSize',7);
    end

    function filterWave(condIdx)
        source = steadyState(condData(condIdx).selectedCond);
        idx = condData(condIdx).waveIdx;
        condData(condIdx).time = condData(condIdx).wavetimes{idx};
        condData(condIdx).wave = condData(condIdx).waveforms{idx};
        condData(condIdx).nt = size(condData(condIdx).wave,2);
        condData(condIdx).iT = min(condData(condIdx).iT,condData(condIdx).nt);
        active = false(size(source.amp));
        filterName = condData(condIdx).filterName;
        if strcmpi(filterName,'none')
            active(:,2:end) = true;
        else
            secondIdx = [];
            if isfield(source,'i1f2'), secondIdx = source.i1f2; end
            filtIdx = determineFilterIndices(filterName,source.freq,source.i1f1,secondIdx);
            active(:,filtIdx) = true;
        end
        if configOptions.applyLowpass
            active(:,source.freq > configOptions.lowpassHz) = false;
        end
        condData(condIdx).activeFreq = active;
        % Unfiltered data always comes from the independently averaged waves.
        % Either a harmonic filter or the low-pass toggle requests filtering.
        if ~strcmpi(filterName,'none') || configOptions.applyLowpass
            filterCfg.activeFreq = active & significantBins(condData(condIdx));
            source.wave = condData(condIdx).wave;
            source.nt = condData(condIdx).nt;
            source.time = condData(condIdx).time;
            condData(condIdx).wave = filterSteadyState(filterCfg,source);
        end
    end
end
