%% Make confusion matrix


bilatRoiList = { 'V1'  'V2D' 'V2V' 'V3A' 'V3D' 'V3V'  'V4'  'LOC' 'MT' ...
    'caudalanteriorcingulate' ...
    'caudalmiddlefrontal' ...
    'frontalpole' ...
    'inferiortemporal' ...
    'isthmuscingulate' ...
    'lateraloccipital' ...
    'lateralorbitofrontal' ...
    'medialorbitofrontal' ...
    'posteriorcingulate' ...
    'rostralanteriorcingulate' ...
    'rostralmiddlefrontal' ...
    'superiorfrontal' ...
    'temporalpole'}


%Ydim{2} = roiList in data
%bilatRoiList = simulated roi's


%Build an index to the corresponding Ydim ROI's for the graphical layout of
%the matrix
for iSim = 1:length(bilatRoiList)
    
    yIdx(iSim) = find(strcmp(bilatRoiList{iSim},Ydim{2}))

end

Ydim{3}

for iSubj = 1:length(Ydim{3}),
    for iSim = 1:length(bilatRoiList)
        
   %create subjectwise confusion matrix, subjxSimulated Source x result
        conMat(iSubj,iSim,:) = abs(squeeze(Yt(1,yIdx,iSubj,4+iSim,1,1,3)))';
   %Normalize by the max source activation     
        conMat(iSubj,iSim,:) = conMat(iSubj,iSim,:)./max(conMat(iSubj,iSim,:));
    end
    
end


%%

imagesc(squeeze(mean(conMat(1,:,:),1)))
set(gca,'ytick',[1:22],'yticklabel',bilatRoiList,'xtick',[1:22])
colorbar 
colormap(gray)
caxis([0 1])



%% Do beamformer conf matrix




%% Set up forward


%mrcDir = '/Volumes/MRI-2/4D2/L1test';
mrcDir = '/Volumes/MRI-1/4D2/RBTX_GV_4D2/GRATING';
mrcDir = '/Volumes/MRI/data/4D2/SEP/mrcProj';
%mrcDir = '/Volumes/MRI/data/4D2/CPT/mrcProj/';
mrcDir = '/Volumes/MRI/data/4D2/JMA_PROJECTS/chopstix/chopstixMRC/';
mrcDir = '/Volumes/MRI/data/4D2/CPT/resolutionSimulations/';


subjList = Ydim_s{3};

tData = zeros(length(subjList),22,22);
lData = tData;
   
    
for iSubj = 1:length(subjList)
    
    
    subjId = subjList{iSubj};

    anatDir = getpref('mrCurrent','AnatomyFolder');
    roiDir = fullfile(anatDir,subjId,'Standard','meshes','ROIs');

    fsDir = getpref('freesurfer','SUBJECTS_DIR');

    fwdFile = fullfile(mrcDir,subjId,'_MNE_',[subjId '-fwd.fif']);
    fwd = mne_read_forward_solution(fwdFile);
    
    % srcFile = fullfile(fsDir,[subjId '_fs4'],'bem',[subjId '_fs4-ico-5p-src.fif']);
    % src = mne_read_source_spaces(srcFile);
    
    src = readDefaultSourceSpace(subjId);
    
    [A Afree] = makeForwardMatrixFromMne(fwd,src);
    
    Anrm = zeros(size(A));
    
    for i=1:length(A);
        
        Anrm(:,i) = A(:,i)./norm(A(:,i));
    end
    
    [funcChunk roiList] = createChunkerFromMeshRoi(roiDir,length(A));

    idx=1;
    roiNames = [];
    for i=1:2:length(roiList),
        roiNames{idx}=roiList(i).name(1:end-6);
        idx=idx+1;
    end

    %Build an index to the corresponding Ydim ROI's for the graphical layout of
    %the matrix
    for iSim = 1:length(bilatRoiList)
        
        myRoiIdx(iSim) = find(strcmp(bilatRoiList{iSim},roiNames));
        
    end
    
    idx = 1;
    funcChunk2 = zeros(size(funcChunk,1),length(myRoiIdx));
    
    for iBilat=myRoiIdx,
        theseRois = 2*(iBilat-1)+1:2*(iBilat-1)+2;
        if ~strcmp(bilatRoiList{idx},roiList(theseRois(1)).name(1:end-6)),
            disp('ERROR')
        end
        
        funcChunk2(:,idx) = sum(funcChunk(:,theseRois),2);
        idx=idx+1;
    end
    
    
    
    for iCond = 1:22,
        
        tD = (squeeze(Y(1,:,iSubj,4+iCond)));
        C=tD'*tD;
    
        [t l tr lr] = lcmvBeamformExtendedSource(A,C+.003*eye(128),funcChunk2);
    
    
        tData(iSubj,iCond,:) = tr./max(tr);
        lData(iSubj,iCond,:) = lr./max(lr);
    end
    
    

end
