%%
% this script maps behavioral data (mat files) onto eye-tracking data (asc
% files) and viceversa, producing a master datasheet ('mat_and_asc_files_AttDeploym.xlsx')
% that contains all files found, with a column that indicates if both kinds 
% of files are present for a given session.

% two other filters are built: one for behavioral data in which sessions
% with less than 15 (out of 20) trials are tagged.
% and a second one for eye data in which session files of less than 7MB are
% tagged.

% the final column of the master datasheet contains an overall filter
% that's TRUE only if all criteria has been passed

%% relevant paths

my_path='/Users/danielrojaslibano/Library/CloudStorage/GoogleDrive-dirl75@gmail.com/.shortcut-targets-by-id/1qQNDY9qvWWjOrFCd88oGEQslFO_Wa08C/LabCode/Magister_NeuroSc/Nunez/AttentionDeployment';
%my_path='G:\.shortcut-targets-by-id\1JD2h14KMMkmsDBGBDB-hW3ZYKKnacMr4\Magister_NeuroSc\Nunez\AttentionDeployment'; %path luzmq


path2saveSheet='/Users/danielrojaslibano/Library/CloudStorage/GoogleDrive-dirl75@gmail.com/My Drive/Matlab_Daniel/Articles/Attentional Deployment/Matlab_scripts_and_functions';
%path2saveSheet='G:\Mi unidad\Tesis\MatLab Tesis\Figuras y resultados'; %path luzma


%if strcmp(computer,'PCWIN64')
    %my_path='G:\.shortcut-targets-by-id\1qQNDY9qvWWjOrFCd88oGEQslFO_Wa08C\LabCode\Magister_NeuroSc\Nunez\AttentionDeployment';

    %path2saveSheet='G:\My Drive\Matlab_Daniel\Articles\Attentional Deployment\Matlab_scripts_and_functions';
%end

% CAR
% my_path='G:\.shortcut-targets-by-id\1qQNDY9qvWWjOrFCd88oGEQslFO_Wa08C\LabCode\Magister_NeuroSc\Nunez\AttentionDeployment';
% path2saveSheet='G:\.shortcut-targets-by-id\1t5cvz3Sd63-V_FV6F2pnX132RtRDWCfy\Attentional Deployment\Matlab_scripts_and_functions';

%%
%two kinds of data files, *.asc and *.mat
dir_asc=dir([my_path filesep '*.asc']);
dir_mat=dir([my_path filesep 'data' filesep '*.mat']);

%get participant's ID from mat files
particip_ID_mat=cell(numel(dir_mat),1); %preallocate
for thisParticipant=1:numel(dir_mat)
    load([dir_mat(thisParticipant).folder filesep dir_mat(thisParticipant).name])
    temp=BehavData.info.Subject_ID;
    %find dashes
    idx_of_dashes=strfind(temp, '-');
    %replace dashes with underscores
    temp(idx_of_dashes)='_';
    particip_ID_mat{thisParticipant}=temp;
end

%get participant's ID from asc files
particip_ID_asc=cell(numel(dir_asc),1); %preallocate
rawFilenames_asc=cell(numel(dir_asc),1); %preallocate
for thisParticipant=1:numel(dir_asc)
    temp=regexp(dir_asc(thisParticipant).name, '.asc', 'split');
    particip_ID_asc{thisParticipant}=temp{1};
    rawFilenames_asc{thisParticipant}=dir_asc(thisParticipant).name;
end


%CHECK FILE SIZE
size_asc_MB=nan(numel(dir_asc),1);
for thisParticipant=1:numel(dir_asc)
    size_asc_MB(thisParticipant)=dir_asc(thisParticipant).bytes/1e6;  
end
%produce tag for too small sessions    
logi_validSessions_asc=size_asc_MB>=6;


% produce tag for incomplete BEHAVIORAL sessions
min_n_of_trials=15; % sessions with less than this n of trials will be tagged
logi_validSessions_behav=true(numel(dir_mat),1); %preallocate
for thisParticipant=1:numel(dir_mat)
    %get behavioral data
    load([dir_mat(thisParticipant).folder filesep dir_mat(thisParticipant).name])
    
    I=BehavData.vars.ResponseIntensity_seq;
    V=BehavData.vars.ResponseValence_seq;
    
    tempI=sum(isnan(I));
    tempV=sum(isnan(V));
    % tag if not enough trials
    if numel(I)- tempV <= min_n_of_trials || numel(I)- tempI<= min_n_of_trials
        logi_validSessions_behav(thisParticipant)=false;
    end
end

% map mat files (behavior) onto asc files (eye-tracking)

%column 1: mat filename
%column 2: asc filename
%column 3: TRUE if map found (ie, a mat & asc pair)
%column 4: TRUE if mat has enough trials
%column 5: TRUE if asc has min size


map_mat2asc=cell(numel(particip_ID_mat),5); %preallocate
for thisParticipant=1:numel(particip_ID_mat)
    
    logi=strcmp(particip_ID_mat{thisParticipant},particip_ID_asc);
    
    if sum(logi==1) % one match, so map
        map_mat2asc{thisParticipant,1}=dir_mat(thisParticipant).name;
        map_mat2asc{thisParticipant,2}=dir_asc(logi).name;
        map_mat2asc{thisParticipant,3}=true;

    elseif sum(logi)==0
        map_mat2asc{thisParticipant,1}=dir_mat(thisParticipant).name;
        map_mat2asc{thisParticipant,2}='Not found';
        map_mat2asc{thisParticipant,3}=false;
    end
    
    map_mat2asc{thisParticipant,4}=logi_validSessions_behav(thisParticipant);
    
end

howMany=size(map_mat2asc,1);
        
%add all asc files that may have been left out
c=0;
for thisParticipant=1:numel(dir_asc)
    
    temp=regexp(dir_asc(thisParticipant).name, '.asc','split');
    temp=temp{1};
    
    logi=strcmp(temp, particip_ID_mat);
    %sum(logi)
    if sum(logi)==0
        c=c+1;
        map_mat2asc{howMany+c,1}='Not Found';
        map_mat2asc{howMany+c,2}=dir_asc(thisParticipant).name;
        map_mat2asc{howMany+c,3}=false;
        map_mat2asc{howMany+c,4}=false;
    end
    
end    

% fill 5th column
for thisParticipant=1:size(map_mat2asc,1)
    
    if strcmp(map_mat2asc{thisParticipant,2}, 'Not found')
        map_mat2asc{thisParticipant,5}=false;
        
    else
        
        logi=strcmp(map_mat2asc{thisParticipant,2},rawFilenames_asc);
        
        if logi_validSessions_asc(logi)
            map_mat2asc{thisParticipant,5}=true;
        else
            map_mat2asc{thisParticipant,5}=false;
        end

    end
    
end

% produce excel output
mat=map_mat2asc(:,1);
asc=map_mat2asc(:,2);
map=cell2mat(map_mat2asc(:,3));
minTrials=cell2mat(map_mat2asc(:,4));
minSizeAsc=cell2mat(map_mat2asc(:,5));

overall=logical(map.*minTrials.*minSizeAsc);

nValidSessions=sum(overall);



MasterTable_Files=table(mat,asc,map,minTrials,minSizeAsc,overall);

%%
writetable(MasterTable_Files, [path2saveSheet filesep 'mat_and_asc_files_AttDeploym.xlsx'])

clearvars -except MasterTable_Files
