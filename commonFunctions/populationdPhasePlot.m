function populationdPhasePlot(phasorData)
cax = gca;

colorList = get(gcf,'DefaultAxesColorOrder');

%colorList = {[0 0 1], [.2 .5 .2],  [1 0 0],[0.01 0.01 0.01],[.5 .2 .5] };

numColorsDefined = length(colorList);

numColors = size(phasorData,2);

colorVec = mod([1:numColors]-1,numColorsDefined)+1;

colorList = colorList(colorVec,:); 

for tmp=1:size(phasorData,2)
    [phi(tmp)] = circ_mean(phasorData(:,tmp)); % mean direction + upper and lower 95% CI
end

% maxComp=max(abs(phasorData));
% maxComp = max(abs(phi));
maxComp = 1;
rmax = max(maxComp,.1);
% rmax = rmax*1.5;
rmax = ceil(10*rmax)/10;
rmax = round(rmax,1,'significant');
rmin = 0;
rticks = 2;
tc = 'k';
ls = '--';
% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% % plot background if necessary
%     if ~ischar(get(cax,'color')),
%        patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
%              'edgecolor',tc,'facecolor',get(cax,'color'),...
%              'handlevisibility','off','parent',cax);
%     end

% draw radial circles
    c82 = cos(78*pi/180);
    s82 = sin(78*pi/180);
    
    rinc = (rmax-rmin)/rticks;
    
    circList = round([(rmin+rinc):rinc:(rmax-rinc)],2,'significant');
    circList = [circList rmax];

    whiteBackH = patch(xunit*circList(end),yunit*circList(end),'w');
    for i=circList;
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off','parent',cax);
        text((i+rinc/30)*c82,(i+rinc/30)*s82, ...
            ['  ' num2str(i,4)],'verticalalignment','bottom',...
            'handlevisibility','off','parent',cax)
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = [0 1.5 3 4.5]*pi/6;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax);

% annotate spokes in degrees
    rt = 1.13*rmax;
    for i = 1:length(th)
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        
        loc = [num2str(th(i)*180/pi +180) '°'];
        text(-rt*cst(i),-rt*snt(i),loc,...
             'horizontalalignment','center',...
             'handlevisibility','off','parent',cax);

        loc = [num2str(th(i)*180/pi) '°'];
        text(rt*cst(i),rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','parent',cax)
    end



%axH = line([-maxComp maxComp; 0 0]',[0 0; -maxComp maxComp]','color','k');
axis tight;
axis equal;
axVal = axis;
%set(axH,'linewidth',1);
tNellip = 30;
tTh = linspace( 0, 2*pi, tNellip )';
SEM = 0;
if SEM
    tNormK = 1 / (length(phasorData)-2);
else % 95% CI
    tNormK = (length(phasorData)-1)/length(phasorData)/(length(phasorData)-2) * finv( 0.95, 2, length(phasorData) - 2 );
end

for iDat = 1:size(phasorData,2)
    
    % mean direction
    r = circ_r(angle(phasorData(:,iDat)));
    phi = circ_mean(angle(phasorData(:,iDat)));
    zm = r*exp(1i*phi);
    centerX = real(zm);
    centerY = imag(zm);
    line([0 centerX], [0, centerY],'Color',colorList(iDat,:), 'linewidth', 2);
    
    % Compute eigen-stuff
    % the idea here is to get two vectors representing the spread of the
    % distribution. 1 vector will be fit to the largest spread and the second
    % will be perpendicular to the 1st vector.
    % tEVec= eigen vector = vector direction (phase) of the 2 vectors. This is
    % only the direction so the second vector will only change its sign (- or
    % +) to be perpendicular to the first one
    % tEVal = eigen values = vector amplitudes. Has only 2 values (one for each
    % vector) with 2 zeros to fit the matrix.
    % tYSubj: one distribution is computed on the column (not row). Multiple
    % clouds should be organised as multiple columns
    [ tEVec, tEVal ] = eig( cov( [ real( phasorData(:,iDat) ), imag( phasorData(:,iDat) ) ] ) );
    % Error/confidence ellipse
    % we start from a circle [ cos(tTh), sin(tTh) ] and ellongate it depending
    % on the spread of the distribution
    tXY = [ cos(tTh), sin(tTh) ] * sqrt( tNormK * tEVal ) * tEVec';
    % the patch has to be drawn at the mean centre
    patch( centerX+tXY(:,1), centerY+tXY(:,2),colorList(iDat,:), 'FaceAlpha', .25, 'EdgeColor', colorList(iDat,:), 'LineWidth', 2 );
    
end

axis(axVal); %plotting the circle can shift the nice circles above.  Reset the value to what it was before the errorbar.

axis off
