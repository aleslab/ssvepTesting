    

load('avTestE2')

for condIdx=1:length(avTest)
    filtIdx = determineFilterIndices( 'low49', avTest(condIdx).freq, avTest(condIdx).i1f1 );
    
    %Create a logical matrix selecting frequency components.
    filtMat = false(size(avTest(condIdx).amp));
    filtMat(:,filtIdx) = true;
    avTest(condIdx).activeFreq = (filtMat)>=1; %Store the filtered coefficients for the spec plot
    cfg.activeFreq =  avTest(condIdx).activeFreq;
    
    [ avTest(condIdx).filteredWave ] = filterSteadyState( cfg, avTest(condIdx) );
end


figure;hold on;
pp=0;
for cc=[5 8]
    pp = pp+1;
    col={'r','b'};
    plot(avTest(cc).time,avTest(cc).filteredWave(23,:),col{pp},'LineWidth',2)
    plot(avTest(cc).time,avTest(cc).wave(23,:),col{pp})
end
legend({'filtered','unfiltered'})