%%% calculates wavelet coherence for CAN regions and physiological signals
%%% includes AR model-based bootstrapping for significance testing
%%% based on Chang and Glover, 2010 "Time–frequency dynamics of resting-state brain connectivity measured with fMRI"

close all
clear all

%% define data parameters, input and output paths
subjects = {'0059','0704','0134','1092','1687','1693','1750','2888',...
    '3046','4901','5262','5571','5931','6173','6293','7002',...
    '7712','7994','8116','8223','8420','8803','8907','9877','9922'};
tasks = {'rest','breathing','coldpressor','rest','breathing'};
sessions = {'02','02','02','03','03'};
NVs = {297,397,777,597,397};
cleanings = {'minimal','wPRF','wCompCor'};
TR = 1.03; %s
Fs = 1/TR;
% input directories
root_dir = "D:\estimation\";
physio_root_dir = "D:\estimation\physio_preproc\";
% output directory
out_root = "D:\CAN_dynamics\WTC\";

%% choose parameters
nboots = 400;
conf_level = 0.95;
plots = false; % choose whether to output figures

%% loop over subjects
for sub_i = 1:length(subjects)
    sub_num = subjects{sub_i};
    for i = 1:length(tasks)
        ses_num = sessions{i};
        task = tasks{i};
        NV = NVs{i};

        %% load in physiological data
        ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task)
        bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);
        path_physio = strcat(physio_root_dir,bids_root,'_physio_and_triggers.mat');                   
        load(path_physio,'time_10','trig_ind_10','HR_10','RF')

        %% extract time
        trig_ind_10 = find(trig_ind_10,NV,"last");            
        t = time_10(trig_ind_10);
        
        for j = 1:length(cleanings)
            cleaning = cleanings{j};

            % prepare output directory
            out_dir = strcat(out_root,cleaning,"/sub-",sub_num,"/ses-",ses_num,"/task-",task,'/');
            [~,~] = mkdir(out_dir);

            %% load in cleaned functional data
            dir_roi_means = strcat(root_dir,'\ROIs\means\',cleaning);
            path_roi_means = strcat(dir_roi_means,'\',ext_root,'_ROI_means.mat'); 
            load(path_roi_means,"roi_means")
    
            %% loop over regions and signals
            n_reg = size(roi_means,2);
            for reg_i = 1:n_reg    
                % pick out CAN region signal
                sig_1 = zscore(roi_means(:,reg_i));
                % fit an AR model (order 1) to the signal
                B_1 = sig_1(1:NV-1)\sig_1(2:NV);
                U_1 = sig_1(2:NV) - sig_1(1:NV-1)*B_1; % residuals
    
                for sig_i = 1:2 %% EXPAND LATER TO INCLUDE EEG?
                    if sig_i == 1
                        sig_phys = "HR";
                        sig_2 = zscore(HR_10(trig_ind_10)).';
                    elseif sig_i == 2
                        sig_phys = "RF";
                        sig_2 = zscore(RF(trig_ind_10));
                    end
    
                    sig_root = strcat('WTC_ROI_',num2str(reg_i),'_x_',sig_phys);
    
                     % fit an AR model (order 1) to the signal
                    B_2 = sig_2(1:NV-1)\sig_2(2:NV);
                    U_2 = sig_2(2:NV) - sig_2(1:NV-1)*B_2; % residuals
    
                    %% calculate wavelet transform coherence
                    [wcoh,wcs,f,coi] = wcoherence(sig_1,sig_2,Fs);
    
                    %% assess significance by bootstrapping
                    wcoh_bs = zeros(size(wcoh,1),NV,nboots);
    
                    % prepare a random array of initializations
                    T0s = randi(NV-1,nboots,1); %%% MAYBE this can be moved outside the loop?
    
                    for nb = 1:nboots
                        T0_b = T0s(nb); % choose the random start value
                    
                        % initialize start values
                        Y_b = zeros(NV,2);
                        Y_b(1,:) = [sig_1(T0_b),sig_2(T0_b)];
                        
                        % prepare randomly-shuffled residuals
                        T_b = randi(NV-1,NV-1,1);
                        U_1_b = U_1(T_b);
                        U_2_b = U_2(T_b);
                    
                        % generate boot-strapped data
                        for t_b = 1:NV-1
                            Y_b(t_b+1,1) = Y_b(t_b,1)*B_1 + U_1_b(t_b);
                            Y_b(t_b+1,2) = Y_b(t_b,2)*B_2 + U_2_b(t_b);
                        end
                    
                        % procede to WTC for boot-strapped time series
                        wcoh_b = wcoherence(Y_b(:,1),Y_b(:,2),Fs);
                        wcoh_bs(:,:,nb) = wcoh_b;
                    end
    
                    % find threshold for statistical significance
                    wcoh_flat = reshape(wcoh_bs,[size(wcoh_bs,1),size(wcoh_bs,2)*size(wcoh_bs,3)]);
                    num_el = round((1-conf_level)*size(wcoh_flat,2)/2);
                    B = maxk(wcoh_flat,num_el,2);
                    conf_thresh = min(B,[],2);% this was around 0.81-0.69 when I used an MAR model
    
                    % create a mask of where original wcoh is significant
                    mask = zeros(size(wcoh));
                    for si = 1:size(wcoh,1)
                        for ti = 1:size(wcoh,2)
                            if wcoh(si,ti) > conf_thresh(si) && f(si) > coi(ti)
                                mask(si,ti) = 1;
                            end
                        end
                    end
    
                    % output plot if desired
                    if plots == true
                        mask_x = repmat(t,size(mask,1),1);
                        mask_y = repmat(f,1,size(mask,2));
                        
                        %mask_surf = pcolor(t,log2(f),mask);
                        %contour_mask = imcontour(mask); % I think that this is returned scaled by axes coordinates
                        %contour_trans = [contour_mask(1,:),]% convert this back to log?
    
                        fig_wtc = figure;                    
                        h = pcolor(t,f,wcoh);
                        h.EdgeColor = "none";
                        ax = gca;
                        ax.YScale='log';
                        yticks = ([0.01 0.05 0.1 0.2 0.3]);
                        ax.YTick = yticks;
                        ax.XLabel.String="Time (s)";
                        ax.YLabel.String="Frequency (Hz)";
                        ax.Title.String = strcat("Wavelet Coherence: CAN",{' '}, ...
                            num2str(reg_i),{' '},'with ',{' '},sig_phys);
                        hcol = colorbar;
                        hcol.Label.String = "Magnitude-Squared Coherence";
                        hold on
                        contour(ax,mask_x,mask_y,mask,1,'white',LineWidth = 1.8);
                        plot(ax,t,coi,"w--",linewidth=2)
                        fig_path = strcat(out_dir,'WTC_plot_',sig_root,'.png');
                        saveas(fig_wtc,fig_path)
    
                        close(fig_wtc)
                    end
    
                    % save results for this pair of signals!
                    out_path = strcat(out_dir,sig_root,'.mat');
                    save(out_path,'wcoh','wcs','f','coi','conf_thresh','mask')
                end
            end
        end
    end
end