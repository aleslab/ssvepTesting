function [newData] = multiplyAxx(data,value)
% multiplication by a given value

% copy all before changing the rest.. 
% so that values that do not change are already copied
newData = data;

% do the addition
newData.wave = data.wave .* value; 
newData.sin = data.sin .* value; 
newData.cos = data.cos .* value;
newData.amp = data.amp .* value;

end
