function [gp] = gpStatT2_DC(sbj,gp)


warning('off','all')
for elec = 2: size(sbj(1,1).sin,1) % elec1 = baseline
for condition = 1: length(sbj)
    for fHz = 2: size(sbj(1,1).sin,2) % fHz is NaN
        clear test
        for sbjNb = 1 : length(sbj(1,:))
            test(sbjNb) = complex(sbj(condition,sbjNb).cos(elec,fHz), sbj(condition,sbjNb).sin(elec,fHz)); % back to complex number
        end
        [pVal, stdDev, confRadius, pT2] = t2circ(complex(test)');
        gp(condition).pval(elec,fHz) = pT2; % because it is group, needs pT2 to allow for covariance instead of pVal 
%         gp(condition).confradius(elec,fHz) = confRadius ; % for the group it is not a circle but an ellipse so should not compute a radius of a circle (as for each individual participant) but the one of an ellipse
%         gp(condition).stderrradius(elec,fHz) = stdDev ; % see above
    end
end
fprintf('elec %d\n',elec)
end

end
        

