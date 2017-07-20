%% harvest data, inverse, & elp-file from mrC project
%clear
%mrCproj = '/Volumes/Denali_4D2/4D2/CPT/groupsDominantEye';
mrCproj = '/Volumes/Denali_4D2/4D2/CPT/groupsAmblyopicEye';

condList = { ...
    'Axx_c001' 'Axx_c002' 'Axx_c003' 'Axx_c004' ...
%    'Axx_c005' 'Axx_c006' 'Axx_c007' 'Axx_c008' ...
    };


groupList = load(fullfile(mrCproj,'SbjGroups.mat'));


 
anatDir = getpref('mrCurrent','AnatomyFolder')
%mrCproj = '/Volumes/MRI-2/4D2/RBTX_GV_4D2/GRATING';condList = {'Axx_c007'};


%invName = 'mneInv_bem_gcv_regu_F1_1_2_3_4_5_6_7_8_9_10_Intermod0.inv';

invName = 'mneInv_bem_nonorm_jma_snr_100.inv'


%condList = {'Axx_c005.mat','Axx_c006.mat'};


subjList = dir([mrCproj '/skeri*']);


toSubj = 'mni0001';

toCtx = readDefaultCortex(toSubj);


clear subjectNameList;

% Create mapping matrices for all subjects, with 3rd order neighbor smoothing
idx = 0;

for iSubj = 1:length(subjList)

    %subjid = 'skeri0047';
    subjid = subjList(iSubj).name

    subjectNameList{iSubj} = subjid;

    mapMtx = makeDefaultCortexMorphMap(subjid,toSubj);
    invFileName = fullfile(mrCproj,subjid,'Inverses',invName);
    fwdFile = fullfile(mrCproj,subjid,'_MNE_',[subjid '-fwd.fif']);
    roiDir = fullfile(anatDir,subjid,'Standard','meshes','ROIs');
    
    %To a second order neighbor smoothing
    fromCtx = readDefaultCortex(subjid);
    fromCtx.uniqueVertices = fromCtx.vertices;
    fromCtx.uniqueFaceIndexList = fromCtx.faces;
    [fromCtx.connectionMatrix] = findConnectionMatrix(fromCtx);
    fromCtx.connectionMatrix = fromCtx.connectionMatrix + speye(length(fromCtx.connectionMatrix)); 
    sumNeighbours=sum(fromCtx.connectionMatrix,2); % Although it should be symmetric, we specify row-summation
    smoothMtx=bsxfun(@rdivide,fromCtx.connectionMatrix,sumNeighbours);
    smoothMtx = smoothMtx*smoothMtx;
    inv = mrC_readEMSEinvFile(invFileName);


	allMapMtx{iSubj} = (mapMtx*smoothMtx)';
    allInv{iSubj}    = inv;
    allCombo{iSubj}  = inv*allMapMtx{iSubj};
%     fwd = mne_read_forward_solution(fwdFile);
%     src = readDefaultSourceSpace(subjid);

%     
 
end
	


%% Load Data
for iSubj = 1:length(subjList)
    
    subjid = subjList(iSubj).name
    
    for iCond = 1:length(condList),
        
        thisCond = condList{iCond};
        E = load(fullfile(mrCproj,subjid,'Exp_MATL_HCN_128_Avg',thisCond));
        
        % allocate space for all data
        if (iCond == 1 && iSubj == 1)
            
            allWave = zeros(length(subjList),length(condList),E.nT,E.nCh);
            Fs= 1000*1/E.dTms;
        end
        
        allWave(iSubj,iCond,:,:) = E.Wave;
        
    end
    
end


%%  Calculate contrasts of interest

%conditions to average together
withinSubjMeanList = { [1 2] 3 4 }
%mean groups to contrast
withinSubjDiffList = { 1 };


allMean = zeros(length(subjList),length(withinSubjMeanList),E.nT,E.nCh);
allDiff = zeros(length(subjList),length(withinSubjDiffList),E.nT,E.nCh);

for iMean = 1:length(withinSubjMeanList),
    
    conds2Mean = withinSubjMeanList{iMean};
    allMean(:,iMean,:,:) = mean(allWave(:,conds2Mean,:,:),2);
end

%         
% for iDiff = 1:length(withinSubjDiffList),
%     
%     conds2Diff = withinSubjDiffList{iDiff};
%     
%     allDiff(:,iDiff,:,:) = allMean(:,conds2Diff(1),:,:)-allMean(:,conds2Diff(2),:,:);
% end
        



%%


baseline = 1:41; 
contrast = [3 ];
thisSubjList = [];

%groupChoice = 'Typicals';
groupChoice = 'Amblyopes';

groupNum  = find(strcmp(groupChoice,{groupList.gGroups.name}))

groupSubjNames = groupList.gGroups(groupNum).members

for iName = 1:length(groupSubjNames),
    
thisSubjList(iName) = find(strcmp(subjectNameList,groupSubjNames{iName}));

end


crossSubCtx = zeros(length(thisSubjList),E.nT,20484);
crossSubSensor = zeros(length(thisSubjList),E.nT,20484);
filtData = zeros(length(thisSubjList),2,E.nT,E.nCh);

subjList

subjIdx = 1;
for iSubj = thisSubjList;%1:length(subjList)

    [num2str(iSubj) ' '  subjectNameList{iSubj}]
    
    thisSub1 = squeeze(allMean(iSubj,contrast(1),:,:));
%     thisSub1 = ftlowpassfilter(thisSub1',Fs,30,15)';
    
    thisSub1 = bsxfun(@minus,thisSub1,mean(thisSub1(baseline,:),1)); %Baseline
    thisSub1 = bsxfun(@minus,thisSub1,mean(thisSub1,2)); %Fix average reference
    filtData(iSubj,1,:,:) = thisSub1;

    
%     thisSub2 = squeeze(allMean(iSubj,contrast(2),:,:));
%     thisSub2 = ftlowpassfilter(thisSub2',Fs,30,15)';
%     
%     thisSub2 = bsxfun(@minus,thisSub2,mean(thisSub2(baseline,:),1)); %Baseline
%     thisSub2 = bsxfun(@minus,thisSub2,mean(thisSub2,2)); %Fix average reference
%     filtData(iSubj,2,:,:) = thisSub2;


%    crossSubCtx(iSubj,:,:) = (thisSub1 - thisSub2)*allCombo{iSubj};
    crossSubCtx(subjIdx,:,:) = thisSub1*allCombo{iSubj};
   subjIdx = subjIdx +1;
   
    
    
end
%

meanCrossSub = squeeze(nanmean(crossSubCtx,1));
%seCrossSub = squeeze(std(crossSubCtx,[],1));
% seCrossSub = seCrossSub./sqrt(size(crossSubCtx,1));
% tCrossSub = meanCrossSub./seCrossSub;


%%
prestimBase = 1:30;
prestimNoise = std(meanCrossSub(prestimBase,:),[],1);
prestimNorm = bsxfun(@rdivide,meanCrossSub,prestimNoise);

%% Plot stuff on the cortex
%thisCol = fromCtx.vertices(:,2);
%toCol = zeros(size( fromCtx.vertices(:,2)));
figure(11);clf
xt = linspace(0,E.dTms*E.nT,E.nT);
interactiveCortexWavePlot(abs(10e6*meanCrossSub),speye(20484),toSubj,xt)




%% Get rois for cross sub average

roiList = { 'LOC-L' 'LOC-R' 'MT-L' 'MT-R' 'V4-L' 'V4-R' 'V1-L' 'V1-R' 'V3A-L' 'V3A-R'};


anatDir = getpref('mrCurrent','AnatomyFolder');
clear allRoi

roiSurfAve = zeros(length(roiList),20484);

for iSubj =1:length(subjList)

    
    subjid = subjList(iSubj).name;
    
    subjRoiDir = fullfile(anatDir,subjid,'Standard','meshes','ROIs');
          
    mapMtx = makeDefaultCortexMorphMap(subjid,toSubj);
    
    for iRoi = 1:length(roiList),
        
        allRoi(iSubj,iRoi) =  load(fullfile(subjRoiDir,roiList{iRoi}));
        thisRoi = zeros(20484,1);
        thisRoi(allRoi(iSubj,iRoi).ROI.meshIndices) = 1;
        roiSurfAve(iRoi,:) = roiSurfAve(iRoi,:) + (mapMtx*thisRoi)';
        
    end
    
    
end


%%

