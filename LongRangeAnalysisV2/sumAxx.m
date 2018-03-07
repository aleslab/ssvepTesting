function [summatedAxx] = sumAxx(Axx1,Axx2)
% just a simple addition

% copy all before changing the rest.. 
% so that values that do not change are already copied
summatedAxx = Axx1;

% do the addition
summatedAxx.wave = Axx1.wave +Axx2.wave;
summatedAxx.sin = Axx1.sin + Axx2.sin;
summatedAxx.cos = Axx1.cos + Axx2.cos;
summatedAxx.amp = sqrt(summatedAxx.sin.^2 + summatedAxx.cos.^2);

end