function [Simult] = reconstructSimult(Axx1,Axx2)
% just a simple addition

% copy all before changing the rest.. 
% so that values that do not change are already copied
Simult = Axx1;

% do the addition
Simult.wave = Axx1.wave +Axx2.wave;
Simult.sin = Axx1.sin + Axx2.sin;
Simult.cos = Axx1.cos + Axx2.cos;
Simult.amp = sqrt(Simult.sin.^2 + Simult.cos.^2);

end