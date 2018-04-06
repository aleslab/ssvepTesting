function [gp] = gpStatT2(sbj,gp)



for elec = 1: size(sbj(1,1).data.sin,1)
for condition = 1: length(sbj(1,:))
    for fHz = 2: size(sbj(1,1).data.sin,2)
        for sbjNb = 1 : length(sbj)
            test(sbjNb) = complex(sbj(sbjNb,condition).data.cos(elec,fHz), sbj(sbjNb,condition).data.sin(elec,fHz)); % back to complex number
        end
        [pVal stdDev confRadius pT2] = t2circ(test');
        gp(condition).pval(elec,fHz) = pT2; % because it is group, needs pT2 to allow for covariance instead of pVal 
%         gp(condition).confradius(elec,fHz) = ; % for the group it is not a circle but an ellipse so should not compute a radius of a circle (as for each individual participant) but the one of an ellipse
%         gp(condition).stderrradius(elec,fHz) = ; % see above
    end
end
fprintf('elec %d\n',elec)
end

end
        