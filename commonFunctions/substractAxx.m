function [subAxx] = substractAxx(Axx1,Axx2)
% just a simple substraction

% copy all before changing the rest.. 
% so that values that do not change are already copied
subAxx = Axx1;

% do the substraction
subAxx.wave = Axx1.wave - Axx2.wave;
subAxx.sin = Axx1.sin - Axx2.sin;
subAxx.cos = Axx1.cos - Axx2.cos;
subAxx.amp = sqrt(subAxx.sin.^2 + subAxx.cos.^2);

end