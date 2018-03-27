function [gp] = gpStatT2(sbj,gp)



for elec = 1: size(sbj(1,1).data.sin,1)
for condition = 1: length(sbj(1,:))
    for fHz = 2: size(sbj(1,1).data.sin,2)
        for sbjNb = 1 : length(sbj)
            test(sbjNb) = complex(sbj(sbjNb,condition).data.cos(elec,fHz), sbj(sbjNb,condition).data.sin(elec,fHz)); % back to complex number
        end
        [pVal stdDev confRadius] = t2circ(test);
        gp(condition).pval(elec,fHz) = pVal;
        gp(condition).confradius(elec,fHz) = confRadius;
        gp(condition).stderrradius(elec,fHz) = stdDev/sqrt(length(test)); % NOT SURE
    end
end
fprintf('elec %d\n',elec)
end

end
        