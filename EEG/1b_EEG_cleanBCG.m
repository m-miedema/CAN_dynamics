% Script to apply ballistocardiogram (BCG) cleaning to EEG data
% Author: Mary Miedema, 2025
% Script should be run inside EEGlab folder; should first run eeglab 
% command to initialize functions
%
% NOTE: it is important to modify the fMRIB script within the eeglab install
% to eliminate the 210 ms delay if the BCG (simulated ECG) channel is used 
% ------------------------------------------------

close all
clear all

subjects = {'0059','0704','0134','1092','1693','1750','2888',...
    '3046','4901','5262','5571','5931','6173','6293','7002',...
    '7712','7994','8116','8223','8420','8803','8907','9877','9922'};
tasks = {'rest','breathing'};

% Empty channels will give a rank warning, turn these off
warning('off','MATLAB:rankDeficientMatrix')

for sub_i = 1:length(subjects)
    subject_id = subjects{sub_i};

    for i = 1:length(tasks)
        task = tasks{i};        

        tic;
        EEG_filename = strcat(subject_id,'_',task,'.set');
        disp(['Now running BCG correction for: ',EEG_filename]);
        EEG = pop_loadset('filename',EEG_filename,'filepath','E:\\eeg_pipeline\\EEGlab_preBCG\\');
        
        EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',40);
        allchan = struct2cell(EEG.chanlocs);
        BCG_ind = find(strcmp(allchan(1,:,:),'BCG'));
        EEG=pop_fmrib_qrsdetect(EEG,BCG_ind,'qrs','no');
        EEG = pop_fmrib_pas(EEG,'qrs','obs',3);
        EEG_filename_out = strcat('E:\\eeg_pipeline\\EEGlab_postBCG\\',subject_id,'_',task,'_bcg_obs_corr.set');
        pop_saveset(EEG, 'filename',EEG_filename_out);
        toc;

    end
end
