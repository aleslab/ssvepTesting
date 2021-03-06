function [] = interactiveTopoSpecPlot(cfg, Axx)
%function [] = interactiveTopoPlot(cfg,Axx)
%
% helpful help here
%



% prepare the layout, this should be done only once
tmpcfg     = removefields(cfg, 'inputfile');
cfg.layout = ft_prepare_layout(tmpcfg);

x = Axx.freq;
data = Axx.Amp;
if ~isfield(Axx,'tCircPval')
    sigFreqs = [];
else
    sigFreqs = Axx.tCircPval<=.01;
end


iElec = 1;
iFr = 1;
%figH= figure; %Render a new figure.

topoAx = subplot(10,1,1:8);
topoH = plotTopo(squeeze(data(:,1)),cfg.layout);

set(gcf,'KeyPressFcn',@keyInput)
colormap(hot);
axis off;
elecVerts = get(topoH,'Vertices');
hold on;
markH = plot(elecVerts(iElec,1),elecVerts(iElec,2),'ko','markersize',10,'linewidth',2);
elecNumH = text(-.05,1.35,num2str(iElec));


freq = x;

%sigFreqs = (handles.data.i1F1:handles.data.i1F1:handles.data.nFr-1)+1;
specAx = subplot(10,1,9:10);
drawSpec();


set(topoH,'ButtonDownFcn',@clickedTopo)
set(topoAx,'ButtonDownFcn',@clickedTopo)




%set(topoAx,'ButtonDownFcn',@specUpdate)

    function drawSpec()
        
        axes(specAx)
        [barH sigH] = pdSpecPlot(freq,data(iElec,:),sigFreqs(iElec,:));
        title(specAx,['Frequency: ' num2str(x(iFr)) ' Hz'])
        
        %Set up the function to call when the plots are clicked on
        set(barH,'ButtonDownFcn',@clickedSpec)
        set(sigH,'ButtonDownFcn',@clickedSpec)
        set(specAx,'ButtonDownFcn',@clickedSpec)
        
    end


    function clickedSpec(varargin)
        
       tCurPoint = get(specAx,'CurrentPoint');
        %set( tHline, 'XData', tCurPoint(1,[1 1]) )
        
        %Get the x location of the lick and find the nearest index
        [distance iFr] = min(abs(x-tCurPoint(1,1)));
        
        if iFr>=1 && iFr<=length(x)

            set(topoH,'facevertexCData',data(:,iFr))
            caxis(topoAx,[0 max(abs(data(:,iFr)))])
            title(specAx,['Frequency: ' num2str(x(iFr)) ' Hz'])

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
            markH = plot(elecVerts(iElec,1),elecVerts(iElec,2),'ko','markersize',10,'linewidth',2);
            set(elecNumH,'String',num2str(iElec));
            title(topoAx,[num2str(iElec) ': ' cfg.layout.label{iElec}])
            drawSpec;

        end
        
    end

    function keyInput(src,evnt)
        
        switch(lower(evnt.Key))
            case 'leftarrow'
                iFr = max(iFr-1,1);
            case 'rightarrow'
                iFr= min(iFr+1,length(x));
                
                
        end
        
        
          set(topoH,'facevertexCData',data(:,iFr))
          caxis(topoAx,[0 max(abs(data(:,iFr)))])
          title(specAx,['Frequency: ' num2str(x(iFr)) ' Hz'])
            
        
    end
        
end
