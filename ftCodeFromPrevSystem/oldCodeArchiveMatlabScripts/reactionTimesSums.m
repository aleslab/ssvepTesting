projectDir = 'X:\data\4D2\CPT\responseLocked';

subjectList = dir(projectDir);

subjectList = subjectList(5:end);

for i=1:length(subjectList)
    RT(i).reaction = load(fullfile(projectDir, subjectList(i).name,'_dev_', 'reactionTimes.mat'));
    
end
    

for i = 1:length(subjectList)
    response(i).name = subjectList(i).name;
    response(i).mean =  mean(RT(i).reaction.reactionTime.buttonLatencySeconds);
    response(i).std = std(RT(i).reaction.reactionTime.buttonLatencySeconds);
    response(i).vect = RT(i).reaction.reactionTime.buttonLatencySeconds;
end

control_ndx = [ 1 2 4 5 6 8 9 12 ];  
ambl_ndx = [ 3 7 10 11 13 14 15 ];

control_mean = 0;
control_std = 0;
control_full_vec = [];
for ndx = 1 : length( control_ndx )
    control_mean = control_mean + 1 / length( control_ndx ) * response(control_ndx(ndx)).mean;
    control_std = control_std + 1 / length( control_ndx ) * response(control_ndx(ndx)).std;
    control_full_vec = [ control_full_vec , response(control_ndx(ndx)).vect ];
end

ambl_mean = 0;
ambl_std = 0;
ambl_full_vec = [];
for ndx = 1 : length( ambl_ndx )
    ambl_mean = ambl_mean + 1 / length( ambl_ndx ) * response(ambl_ndx(ndx)).mean;
    ambl_std = ambl_std + 1 / length( ambl_ndx ) * response(ambl_ndx(ndx)).std;
    ambl_full_vec = [ ambl_full_vec , response(ambl_ndx(ndx)).vect ];
end

ambl_eye_Dir = 'X:\data\4D2\CPT\amplyopes';

ambl_subjList = dir(ambl_eye_Dir);

ambl_subjList = ambl_subjList(4:end);

for i=1:length(ambl_subjList)
    ambl_RT(i).reaction = load(fullfile(ambl_eye_Dir, ambl_subjList(i).name,'_dev_', 'reactionTimes.mat'));
    
end

for i = 1:length(ambl_subjList)
    ambl_response(i).name = ambl_subjList(i).name;
    ambl_response(i).mean =  mean(ambl_RT(i).reaction.reactionTime.buttonLatencySeconds);
    ambl_response(i).std = std(ambl_RT(i).reaction.reactionTime.buttonLatencySeconds);
    ambl_response(i).vect = ambl_RT(i).reaction.reactionTime.buttonLatencySeconds;
end

ambl_eye_mean = 0;
ambl_eye_std = 0;
ambl_eye_full_vec = [];
for idx = 1 : length( ambl_subjList )
    ambl_eye_mean = ambl_eye_mean + 1 / length( ambl_subjList ) * ambl_response(idx).mean;
    ambl_eye_std = ambl_eye_std + 1 / length( ambl_subjList ) * ambl_response(idx).std;
    ambl_eye_full_vec = [ ambl_eye_full_vec , ambl_response(idx).vect ];
end


% for i = 1:length(subjectList)
%     respMat = zeros(length(subjectList),2);
%     respMat(i,1) =  response{i}.mean;
%     respMat(i,2) = response{i}.std;
%     7
% end


% Jacob, you can use these lines to display the bar graph corresponding to
% the averaged time responses (the small bars are the standard errors).
figure
bar( [control_mean ; ambl_mean ; ambl_eye_mean ] )
hold on
bar( [control_mean ; 0 ] )
bar( [0 ; ambl_mean ; 0 ]);
errorbar([1 2 3]',[control_mean ; ambl_mean ; ambl_eye_mean ],[control_std ; ambl_std ; ambl_eye_std ] ./ std([ length( control_ndx ) length( ambl_ndx ) length(ambl_subjList) ]') )
% errorbar([1 2]',[control_mean ; ambl_mean ],[control_std ; ambl_std ] ./ std([ length( control_ndx ) length( ambl_ndx ) ]') )


