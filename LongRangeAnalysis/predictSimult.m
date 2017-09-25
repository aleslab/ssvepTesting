function [Simult] = predictSimult(dataAxx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULTANEOUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 is the original long range
% 2 the prediction long range
% 3 is the original short range
% 4 is the prediction short range
% 5 the difference long range (original - prediction) -> take the absolute
% value for the amp otherwise won't plot the negative 
% 6 the difference short range

Simult(1) = dataAxx(4);
Simult(3) = dataAxx(8);

% stuff needed but does not change
Simult(2).time =  dataAxx(1).time ;
Simult(2).nt =  dataAxx(1).nt;
Simult(2).freq =  dataAxx(1).freq; % in this case, same frequencies
% have to find out later, for now just put NaN
Simult(2).pval =  NaN(size(dataAxx(1).pval));
Simult(2).confradius =  NaN(size(dataAxx(1).confradius));
% do the addition
Simult(2).wave = dataAxx(2).wave + dataAxx(3).wave;
Simult(2).sin = dataAxx(2).sin + dataAxx(3).sin;
Simult(2).cos = dataAxx(2).cos + dataAxx(3).cos;
Simult(2).amp = sqrt(Simult(2).sin.^2 + Simult(2).cos.^2);


% do the addition
Simult(4).wave = dataAxx(6).wave + dataAxx(7).wave;
Simult(4).sin = dataAxx(6).sin + dataAxx(7).sin;
Simult(4).cos = dataAxx(6).cos + dataAxx(7).cos;
Simult(4).amp = sqrt(Simult(4).sin.^2 + Simult(4).cos.^2);


% do the difference
Simult(5).wave = Simult(1).wave-Simult(2).wave;
Simult(6).wave = Simult(3).wave-Simult(4).wave;
Simult(5).sin = Simult(1).sin-Simult(2).sin;
Simult(6).sin = Simult(3).sin-Simult(4).sin;
Simult(5).cos = Simult(1).cos-Simult(2).cos;
Simult(6).cos = Simult(3).cos-Simult(4).cos;
Simult(5).amp = sqrt(Simult(5).sin.^2 + Simult(5).cos.^2);
Simult(6).amp = sqrt(Simult(6).sin.^2 + Simult(6).cos.^2);

% stuff needed for plotting but does not change
for nb=1:6
    Simult(nb).time =  dataAxx(1).time ;
    Simult(nb).nt =  dataAxx(1).nt;
    Simult(nb).freq =  dataAxx(1).freq; % in this case, same frequencies
    Simult(nb).pval =  NaN(size(dataAxx(1).pval));
    Simult(nb).confradius =  NaN(size(dataAxx(1).confradius));
end


%%% add labels
% 1 is the original long range
Simult(1).label = 'original LR';
% 2 the prediction long range
Simult(2).label = 'prediction LR';
% 3 is the original short range
Simult(3).label = 'original SR';
% 4 is the prediction short range
Simult(4).label = 'prediction SR';
% 5 the difference long range (original - prediction) 
Simult(5).label = 'difference LR';
% 6 the difference short range
Simult(6).label = 'difference SR';

end