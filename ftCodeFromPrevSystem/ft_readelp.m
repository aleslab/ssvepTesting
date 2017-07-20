%%readelp
elecfile = 'GD_AmblyCPT4_11302009.elp';
elp = readelp(elecfile);

temppnt = [elp.X; elp.Y; elp.Z]';
elec.pnt(1:128,:) = temppnt(5:end,:);
elec.pnt(129,:) = temppnt(4,:);

templabel = {elp.labels};
elec.label = templabel(5:end);
elec.label{129} = templabel{4};

cfg1.elec = elec;


    