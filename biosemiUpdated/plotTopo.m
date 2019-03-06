function varargout = plotTopo(data,layout)
%plotTopo - Plots data on arbitrary layout
%
%

if ~isvector(data)
    error('Data input must be a vector')
end

data = squeeze(data); %Remove any singleton dimensions
data = data(:);  %Make column matrix
dataSz = size(data); %Fin

nChan = dataSz(1); %Number of channels in data 

if ischar(layout)
    tmpcfg.layout = layout
    layout = ft_prepare_layout(tmpcfg);
end


if isstruct(layout) && isfield(layout,'cfg') %If cfg field is present layout is FT format
       
    tEpos = [ layout.pos(1:nChan,1), layout.pos(1:nChan,2), zeros(nChan,1) ];
    
    tFaces = delaunay(tEpos(1:nChan,1),tEpos(1:nChan,2));

end


patchList = findobj(gca,'type','patch');
netList   = findobj(patchList,'UserData','plotTopo');


if isempty(netList),    
    handle = patch( 'Vertices', [ tEpos(1:nChan,1:2), zeros(nChan,1) ], ...
        'Faces', tFaces,'EdgeColor', [ 0.5 0.5 0.5 ], ...
        'FaceColor', 'interp');
    axis equal;
    axis off;
    curAx= axis;
    axis(curAx*1.1); % increases axis limit by 10%
else
    handle = netList;
end

set(handle,'facevertexCdata',data,'linewidth',1,'markersize',20,'marker','.');
set(handle,'userdata','plotTopo');

% colormap(jmaColors('usadarkblue'));

if nargout >= 1
varargout{1} = handle;
end
