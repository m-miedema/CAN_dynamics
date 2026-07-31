% sort FD, vigilance, HR, and RF values by CAP states

% choose script options
z_scoring = 0; % toggles whether vigilance, fd, and hr are z-scored within subjects
make_overview_fig = 0;
make_comparison_fig = 1;
ica_choice = 1;
states_to_explore = {[2:15],[2:15]};
vig_dir = "E:\eeg_pipeline\EEG_measures\";
mot_root = "E:\estimation\";
physio_root_dir = "E:\estimation\physio_preproc\";
dFC_dir = "E:\CAN_dynamics\dFC_CAP\";
fig_dir = "E:\CAN_dynamics\CAP_states\";
tasks = {'rest','breathing'};
dFC_lens = {597,397};
vig_lens = {595,395};
sessions = {'03','03'};
comp_types = {'vigilance','fd','rf','hr'};
cleanings = {'minimal','wPRF','wCompCor'};
ica_crit = {'aggressive','conservative'};
n_state_comps = 14; % number of states that are compared

TR = 1.03;

pipeline_colors = [[0.3 0.2 0.4];[0.356 0.456 0.758];[0.36 0.6 0.25]];
pipeline_labels = {'M+RETROICOR','+PRF','+CompCor'};


%% create a place to store information about coincident differences
diffs_fd = {};
diffs_hr = {};
diffs_rf = {};
diffs_vig = {};

% also store differences between centroids associated with variable changes
centroid_diffs_fd = {{},{},{},{},{},{}};
centroid_diffs_hr = {{},{},{},{},{},{}};
centroid_diffs_rf = {{},{},{},{},{},{}};
centroid_diffs_vig = {{},{},{},{},{},{}};
centroid_diff_avgs = zeros(numel(comp_types),2,numel(cleanings),n_state_comps,31);

%% function for z-scoring which can handle NaNs
myzscore = @(x) (x - nanmean(x)) ./ nanstd(x);

%% loop over comparison variables
for m = 1%:numel(comp_types)
    comp_type = comp_types{m};
    % p_values_all = {};
    % n_states_all = {};
    %% loop over scans
    for i = 1:numel(tasks)
        ses_num = sessions{i};
        task = tasks{i};
        vig_len = vig_lens{i};
        dFC_len = dFC_lens{i};
        exp_states = states_to_explore{i};
    
        if strcmp('rest',task)
            subjects = {'0059','0704','0134','1092','1693','1750','2888',...
        '3046','5262','5571','5931','6173','6293','7002',...
        '7712','7994','8116','8223','8420','8803','8907','9877','9922'};
        elseif strcmp('breathing',task)
            subjects = {'0059','0704','1092','1693','1750','2888',...
         '3046','4901','5262','5571','5931','6173','6293','7002',...
        '7712','7994','8116','8420','8907','9877'};
        end
        p_values_all = {};
        n_states_all = {};
        %% loop over pipelines
        for c_i = 1:numel(cleanings)
            cleaning=cleanings{c_i};   
    
            % create arrays to hold vigilance and state labels
            dFC_labels = zeros(numel(subjects),n_state_comps,dFC_len);
            vigil_ts = zeros(numel(subjects),numel(ica_crit),vig_len);
            hr_ts = zeros(numel(subjects),dFC_len);
            rf_ts = zeros(numel(subjects),dFC_len);
            fd_ts = zeros(numel(subjects),dFC_len-1);
            % create arrays to hold clustering statistics
            dwell_t_cell = {};
            transition_p_cell = {};
            for n_i = 2:(1+n_state_comps)
                dwell_t_cell{end+1} = zeros(numel(subjects),n_i);
                transition_p_cell{end+1} = zeros(numel(subjects),n_i,n_i);
            end
            % array to track how many states were specified
            n_states = [];
    
            %% loop over subjects
            for s_i = 1:numel(subjects)
                sub_num = subjects{s_i};
    
                %% loop over number of CAP states used in dFC model
                for n_i = 0:n_state_comps-1
                    % load in data for this subject
                    dFC_path = strcat(dFC_dir,cleaning,'\sub-',sub_num,'_ses-',ses_num,'_task-',task,'_n-',int2str(n_i),'_CAP_labels.mat');
                    load(dFC_path)
                    % store time series
                    dFC_labels(s_i,n_i+1,:) = FCS_idx_array;
                    % get the number of states actually specified
                    n_states_i = size(FCS_proba,2);
    
                    %% calculate transition probability
                    [FCS_no_repeats, ~, n_level] = unique(FCS_idx_array);
                    FCS_idx_array = FCS_idx_array + 1; % switch from Python to Matlab indexing
                    counts = accumarray([FCS_idx_array(1:end-1)', FCS_idx_array(2:end)'], 1, [n_states_i, n_states_i]);
                    transition_p = counts ./ sum(counts, 2);
                    transition_p(isnan(transition_p)) = 0; 
                    % store in cell array
                    transition_p_cell{n_states_i-1}(s_i,:,:) = transition_p;
    
                    %% calculate dwell time
                    change_states = find(diff(FCS_idx_array) ~= 0);
                    change_state_pts = [0, change_states, length(FCS_idx_array)];
                    dwell_t = diff(change_state_pts)*TR;
                    % find states at start of transition
                    states = FCS_idx_array(change_state_pts(1:end-1) + 1);
                    % calculate the mean dwell time for each state
                    mean_dwell_t = zeros([n_states_i 1]);
                    for state = 1:n_states_i
                        state_durations = dwell_t(states == state);
                        mean_dwell_t(state) = mean(state_durations);
                    end
                    % store in cell array
                    dwell_t_cell{n_states_i-1}(s_i,:) = mean_dwell_t;
    
                    % store number of states actually specified
                    if s_i == 1
                        n_states(end+1) = n_states_i;
                    end
                end
    
                if strcmp('vigilance',comp_type)
                    %% loop over ICA approach and load in vigilance
                    for ica_i = 1:numel(ica_crit)
                        ica = ica_crit{ica_i};
                        vig_path = strcat(vig_dir,'sub-',sub_num,'_task-',task,'_ICA-',ica,'_vigilance.mat');
                        load(vig_path)
                        vigilance_good_epochs = vigilance;
                        vigilance_good_epochs(bad_epochs+1) = NaN; % convert from Python index
                        if z_scoring == 1
                            Z = myzscore(vigilance_good_epochs');
                            vigilance_good_epochs = Z';
                        end
                        vigil_ts(s_i,ica_i,:) = vigilance_good_epochs;
                    end
                elseif strcmp('fd',comp_type)
                    % load in head motion and calculate framewise displacement
                    mot_path = strcat(mot_root,'\motion\',"sub-",sub_num,"_ses-",ses_num,"_",task,'_mot.par');
                    mot_params = load(mot_path);
                    mot_d = abs(diff(mot_params));
                    if z_scoring == 0
                        fd_ts(s_i,:) = sum(mot_d(:,4:6),2) + 50*(pi/180)*sum(mot_d(:,1:3),2);
                    elseif z_scoring == 1
                        fd_ts(s_i,:) = zscore(sum(mot_d(:,4:6),2) + 50*(pi/180)*sum(mot_d(:,1:3),2));
                    end
                else
                    % load in physiological signals
                    ext_root = strcat("sub-",sub_num,"_ses-",ses_num,"_task-",task);
                    bids_root = strcat("sub-",sub_num,"/ses-",ses_num,"/func/",ext_root);
                    path_physio = strcat(physio_root_dir,bids_root,'_physio_and_triggers.mat');                   
                    load(path_physio,'trig_ind_10','HR_10','RF')
                    trig_ind_10 = find(trig_ind_10,dFC_len,"last");
                    rf_ts(s_i,:) = zscore(RF(trig_ind_10)).';                    
                    if z_scoring == 0
                        hr_ts(s_i,:) = HR_10(trig_ind_10).';
                    elseif z_scoring == 1
                        hr_ts(s_i,:) = zscore(HR_10(trig_ind_10)).';
                    end
                end
        
            end
    
            p_values = [];
    
            %% loop over number of CAP states
            for n_i = 1:n_state_comps
    
                n_act = n_states(n_i);
    
                dFC_labels_i = squeeze(dFC_labels(:,n_i,:)); %switch from python to Matlab indexing
    
                % now get comparatory time series
                if strcmp('vigilance',comp_type)
                    comp_i = squeeze(vigil_ts(:,ica_choice,:));
                    dFC_labels_i = dFC_labels_i(:,3:end);
                elseif strcmp('fd',comp_type)
                    comp_i = fd_ts;
                    dFC_labels_i = dFC_labels_i(:,2:end);
                elseif strcmp('rf',comp_type)
                    comp_i = rf_ts;
                elseif strcmp('hr',comp_type)
                    comp_i = hr_ts;
                end
    
                % run a one-way anova
                [p,tbl,stats] = anova1(comp_i(:),dFC_labels_i(:),"off");
                %[c,m] = multcompare(stats,'Display', 'off'); % last column of c is the p-value for comparison
                p_values(end+1) = p;
    
                % calculate the % of time spent in each state for each subject
                frac_i = zeros(numel(subjects),n_act);
                for state_i = 0:(n_act-1)
                    counts = zeros(size(dFC_labels_i));
                    counts(dFC_labels_i==state_i)=1;
                    frac_i(:,state_i+1) = sum(counts,2)/vig_len;
                end
                mean_percent_i = mean(frac_i);
                xs = repmat(1:n_act,numel(subjects),1);
    
                CAP_states_path = strcat(dFC_dir,cleaning,'\\FCS_states_ses-',ses_num,'_task-',task,'_n-',int2str(n_act),'.mat');
                load(CAP_states_path)
    
                %% plot an overview for all states
    
                if make_overview_fig == 1
        
                    overview_fig = figure(1);
                    n_layout = ceil(n_act/5);
                    set(overview_fig, 'Position', [200, 100, 1000, 200+120*n_layout]);            
                    outer_t = tiledlayout(1+n_layout,1);%,"TileSpacing","tight");
                    inner_t = tiledlayout(outer_t,n_layout,5);
                    inner_t.Layout.Tile = 1;
                    inner_t.Layout.TileSpan = [n_layout 1]; 
                    for state_i = 1:n_act
                        nexttile(inner_t);
                        imagesc(squeeze(FCS_states(state_i,:,:)))
                        axis square;
                        xticks([])
                        yticks([])
                        title(int2str(state_i));
                        c_max = max(abs(FCS_states(state_i,:,:)),[],'all'); 
                        clim([-c_max c_max]);
                        colorbar
                    end
    
                    inner_t_2 = tiledlayout(outer_t,1,3);
                    inner_t_2.Layout.Tile = n_layout+1;
                    % add clustering statistics
                    ax1 = nexttile(inner_t_2);
                    boxplot(frac_i(:),xs(:),'PlotStyle','compact')
                    ylabel({'Fractional'; 'occupancy'})
                    dwell_i = dwell_t_cell{n_act-1};
                    ax2 = nexttile(inner_t_2);
                    boxplot(dwell_i(:),xs(:),'PlotStyle','compact')
                    ylabel("Dwell time (s)")
                    nexttile(inner_t_2)
                    trans_p_i_all = transition_p_cell{n_act-1};
                    trans_p_i = squeeze(mean(trans_p_i_all,"omitmissing"));
                    imagesc(trans_p_i)
                    axis square
                    ylabel({'Transition'; 'probability'})
                    colorbar
    
                    fig_path = strcat(fig_dir,cleaning,'_',task,'_',int2str(n_act),'.png');
                    saveas(overview_fig,fig_path)
                end
        
            end
    
            p_values_all{end+1} = p_values;
            n_states_all{end+1} = n_states;
    
            %% now do specific explorations of state x (electro)physiology
            if make_comparison_fig == 1
                for act_i = 1:numel(exp_states)
        
                    n_act = exp_states(act_i);
                    n_i = find(n_states == n_act);
        
                    % re-load CAP state labels
                    dFC_labels_i = squeeze(dFC_labels(:,n_i,:));
                    if strcmp('vigilance',comp_type)
                        dFC_labels_i = dFC_labels_i(:,3:end);
                    elseif strcmp('fd',comp_type)
                        dFC_labels_i = dFC_labels_i(:,2:end);
                    end

                    % re-load CAP states
                    CAP_states_path = strcat(dFC_dir,cleaning,'\\FCS_states_ses-',ses_num,'_task-',task,'_n-',int2str(n_act),'.mat');
                    load(CAP_states_path)
        
                    % create a table with the median values
                    medians = zeros([n_act 1]);
                    for med_i = 1:n_act
                        medians(med_i) = median(comp_i(dFC_labels_i(:)+1==med_i),'omitnan');
                    end
                    
                    % re-run one-way anova
                    [p,tbl,stats] = anova1(comp_i(:),dFC_labels_i(:),"off");
                    [c,means] = multcompare(stats,'Display', 'off'); % last column of c is the p-value for comparison
        
                    % first find significant differences
                    sig_thresh = 0.05/size(c,1);
                    sig_diffs = find(c(:,6)<sig_thresh); % find rows that are significantly different

                    

                    centroid_diffs = zeros(numel(sig_diffs),size(centroids,2));
                    state_diffs = zeros(1,numel(sig_diffs));
                    sig_diff_labels = {};
                    
                    % prepare to track differences
                    diffs_i = zeros(n_act);

                    % prepare to deal with states that don't contain any samples
                    empty_states = [];

                    for pair_i = 1:numel(sig_diffs)
                        sig_diff = sig_diffs(pair_i);
                        % find the states that are different
                        state_a = c(sig_diff,1);
                        state_b = c(sig_diff,2);
                        % find the medians
                        med_a = medians(state_a);
                        med_b = medians(state_b);
                        if isnan(med_a) || isnan(med_b)
                            % one of these states has no samples, don't
                            % include it
                            empty_states = [empty_states,pair_i];
                            continue
                        end
                        if med_a > med_b
                            centroid_diffs(pair_i,:) = centroids(state_a,:)-centroids(state_b,:);
                            sig_diff_labels{end+1} = strcat(num2str(state_a),'-',num2str(state_b));
                            diffs_i(state_a,state_b) = 1;
                            diffs_i(state_b,state_a) = -1;
                        elseif med_b > med_a
                            centroid_diffs(pair_i,:) = centroids(state_b,:)-centroids(state_a,:);
                            sig_diff_labels{end+1} = strcat(num2str(state_b),'-',num2str(state_a));
                            diffs_i(state_a,state_b) = -1;
                            diffs_i(state_b,state_a) = 1;
                        end
                        state_diffs(pair_i) = abs(med_a - med_b);                        
                    end  

                    % deal with states that don't contain any samples
                    centroid_diffs(empty_states,:) = [];
                    state_diffs(empty_states) = [];

                    % append results to the appropriate main cell array
                    if strcmp('vigilance',comp_type)
                        diffs_vig{end+1} = diffs_i;
                        centroid_diffs_vig{c_i+3*(i-1)} = [centroid_diffs_vig{c_i+3*(i-1)};centroid_diffs];
                    elseif strcmp('fd',comp_type)
                        diffs_fd{end+1} = diffs_i;
                        centroid_diffs_fd{c_i+3*(i-1)} = [centroid_diffs_fd{c_i+3*(i-1)};centroid_diffs];
                    elseif strcmp('hr',comp_type)
                        diffs_hr{end+1} = diffs_i;
                        centroid_diffs_hr{c_i+3*(i-1)} = [centroid_diffs_hr{c_i+3*(i-1)};centroid_diffs];
                    elseif strcmp('rf',comp_type)
                        diffs_rf{end+1} = diffs_i;
                        centroid_diffs_rf{c_i+3*(i-1)} = [centroid_diffs_rf{c_i+3*(i-1)};centroid_diffs];
                    end                    

                    % create a weighted average
                    mean_centroid_diff = state_diffs*centroid_diffs/sum(state_diffs);
                    sig_diff_labels{end+1} = 'avg';
                    centroid_diff_avgs(m,i,c_i,n_act-1,:) = mean_centroid_diff;               

                    % plot the difference between centroids
                    centroid_diff_fig = figure;
                    imagesc([centroid_diffs;mean_centroid_diff])
                    yticks(1:numel(sig_diffs)+1)
                    yticklabels(sig_diff_labels)
                    colorbar
                    centroid_diff_title = strcat(num2str(i),'-',cleaning,":",{' '},num2str(n_act),', ',comp_type);
                    title(centroid_diff_title)
                    if z_scoring == 0
                        fig_path = strcat(fig_dir,cleaning,'_',task,'_',int2str(n_act),'_',comp_type,'_centroid_diffs.png');                    
                    elseif z_scoring == 1
                        fig_path = strcat(fig_dir,cleaning,'_',task,'_',int2str(n_act),'_',comp_type,'_centroid_diffs_Z.png');                    
                    end
                    saveas(centroid_diff_fig,fig_path)
                    close(centroid_diff_fig)
                    
                    comparison_fig = figure;
                    title(cleaning)          
                    t = tiledlayout(1,3,"TileSpacing","tight");
                    nexttile([1 2])
                    boxplot(comp_i(:),dFC_labels_i(:)+1,'PlotStyle','compact');
                    xlabel('CAP state')
                    ylabel(comp_type)
                    hold on
        
                    y_max = max(comp_i(:));
                    y_min = min(comp_i(:));
                    if ~isempty(sig_diffs)
                        axis([xlim    y_min  ceil(y_max*(1.02 + 0.1*numel(sig_diffs)))])
                    end
        
                    for bar_i = 1:numel(sig_diffs)
                        sig_diff = sig_diffs(bar_i);
                        % find the states that are different
                        state_a = c(sig_diff,1);
                        state_b = c(sig_diff,2);
                        y_level = 1 + 0.06*bar_i;
                        plot([state_a, state_b], [1 1]*y_max*y_level, '-k',  mean([state_a, state_b]), max(y_max)*(y_level+0.02), '*k')
                    end
                    hold off
    
                    
                    
                    T = table((1:n_act)', medians, 'VariableNames', {'State', 'Median'});
    
                    % Get the table in string form.
                    TString = evalc('disp(T)');
                    % Use TeX Markup for bold formatting and underscores.
                    TString = strrep(TString,'<strong>','\bf');
                    TString = strrep(TString,'</strong>','\rm');
                    TString = strrep(TString,'_','\_');
                    % Get a fixed-width font.
                    FixedWidth = get(0,'FixedWidthFontName');
                    % Output the table using the annotation command.
                    ax = nexttile;
                    axis off
                    ah = annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
                        'FontName',FixedWidth,'Units','Normalized','Position',[0.6 -0.3 1 1]);
                    ah.LineStyle = 'none';
                    if z_scoring == 0
                        fig_path = strcat(fig_dir,cleaning,'_',task,'_',int2str(n_act),'_',comp_type,'_comparison.png');                    
                    elseif z_scoring == 1
                        fig_path = strcat(fig_dir,cleaning,'_',task,'_',int2str(n_act),'_',comp_type,'_comparison_Z.png');
                    end
                    saveas(comparison_fig,fig_path)
                    close(comparison_fig)
        
                end
            end
    
    
        end
    
        %% create a table to compare
        for c_i = 1:numel(cleanings)
            T = table(n_states_all{c_i}',p_values_all{c_i}', 'VariableNames', {'n_states',cleanings{c_i}});
            if c_i == 1
                masterTable = T;
            else
                masterTable = outerjoin(T,masterTable,'MergeKeys', true);
            end
        end
        if z_scoring == 0
            table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,".xlsx");
        elseif z_scoring == 1
            table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,"_Z.xlsx");
        end
        writetable(masterTable, table_path)
    end
end


%% make a figure of all p-values
p_overview_fig = figure('Position',[360.0000    2.3333  884  638.6667]);
task_symbol = {'o','^'};
if z_scoring == 0
    comp_type_labels = {'Vigilance','Framewise displacement','Respiratory flow (z-scored)','Heart rate'};
elseif z_scoring == 1
    comp_type_labels = {'Vigilance (z-scored)','Framewise displacement (z-scored)','Respiratory flow (z-scored)','Heart rate (z-scored)'};
end
t = tiledlayout(2,2,'TileSpacing','Compact');
title(t,'Significance of differences in variable distributions between CAP states')
for m = 1:numel(comp_types)
    comp_type = comp_types{m};
    comp_type_label = comp_type_labels{m};
    nexttile
    yline(-log10(0.001), ':r','LineWidth',1.3); % significance threshold
    yline(-log10(0.01), '--r','LineWidth',1.3); % significance threshold
    yline(-log10(0.999), 'white','LineWidth',0.1); % dummy variable
    title(comp_type_label)
    xlabel('Number of states','FontSize',10)
    ylabel('-log10(p)','FontSize',10)
    hold on
    for i = 1:numel(tasks)
        task = tasks{i};
        if z_scoring == 0
            table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,".xlsx");
        elseif z_scoring == 1
            table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,"_Z.xlsx");
        end
        % if m == 1
        %     %table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,"_HRnoZ.xlsx"); % NOT ZSCORING HR
        %     table_path = strcat("CAP_sorted_p_values_",task,"_",comp_type,"_vigZ.xlsx")
        % end
        pTable = readtable(table_path);

        % plot p-values
        n_states = pTable.n_states;
        %pipeline_colors = [[0.3 0.1 0.4];[0.356 0.456 0.758];[0.3 0.6 0.25]];
        for c_i = 1:numel(cleanings)
            cleaning=cleanings{c_i};
            logP = -log10(pTable.(cleaning));
            plot(n_states,logP,task_symbol{i},'MarkerEdgeColor',pipeline_colors(c_i,:),'LineWidth',1.6)
        end
    end
    xlim([1 16])
end
legend_labels = {'p < 0.001','p < 0.01','Rest, minimally denoised','Rest, denoised with PRF','Rest, denoised with CompCor','Breathing, minimally denoised','Breathing, denoised with PRF','Breathing, denoised with CompCor'};
legend_labels = {'p < 0.001','p < 0.01',' ','Rest, M+RETROICOR','Rest, +PRF','Rest, +CompCor','Breathing, M+RETROICOR','Breathing, +PRF','Breathing, +CompCor'};
legend(legend_labels,'Location','southoutside','NumColumns',3)
if z_scoring == 0
    fig_path = strcat(fig_dir,'all_p_vals.png');
elseif z_scoring == 1
    fig_path = strcat(fig_dir,'all_p_vals_z.png');
end
saveas(p_overview_fig,fig_path)

%% coincident variable differences
% make colormap
n = 50;
cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
cmap = flip([cmap1; cmap2(2:end, :)]);

coincident_fig = figure('Position',1.0e+03 *[0.0063    0.1957    1.2727    0.3007]);
t = tiledlayout(3,14,"TileSpacing","compact");
title(t,'Coincidence of significant differences in physiological variables during CAP states')
diff_mats = zeros(4,4,3,n_state_comps);
for c_i = 1:3
    for n_i = 1:n_state_comps 
        n_act = n_states(n_i);
        diff_mat = nan(4);
        for i = 1:2
            if i == 1
                i_ind = find(triu(ones(4),1));
            elseif i == 2
                i_ind = find(tril(ones(4),-1));
            end
            ind = n_i + (c_i-1)*n_state_comps + (i-1)*(3*n_state_comps);
            % create a coincidence matrix
            diff_fd = diffs_fd{ind};
            diff_hr = diffs_hr{ind};
            diff_rf = diffs_rf{ind};
            diff_vig = diffs_vig{ind};
            tri_ind = find(triu(ones(size(diff_vig)),1));
            diff_vec = [diff_vig(tri_ind),diff_fd(tri_ind),diff_rf(tri_ind),diff_hr(tri_ind)];
            if size(diff_vec,1) == 1
                diff_mat_i = diff_vec'*diff_vec; % outer product
            else
                diff_mat_i = corrcoef(diff_vec);
            end
            diff_mat(i_ind) = diff_mat_i(i_ind);
            
        end
        diff_mats(:,:,c_i,n_i) = diff_mat;
        nexttile
        h_im = imagesc(diff_mat,[-1 1]);
        nan_spots = ~isnan(diff_mat);
        set(gca, 'Color', [0.5 0.5 0.5]); % set background to gray
        set(h_im, 'AlphaData', nan_spots);
        colormap(cmap)
        yticks([])
        xticks([])
        if c_i == 1
            tile_title = strcat(num2str(n_act),{' states'});
            title(tile_title)
        end
        if n_act == 2
            ylabel_txt = {pipeline_labels{c_i},' '};
            ylabel(ylabel_txt)
        end

    end
end
cb = colorbar;
ylabel(cb,'Coincidence')
if z_scoring == 0
    fig_path = strcat(fig_dir,'CAP_states_coincidence.png');
elseif z_scoring == 1
    fig_path = strcat(fig_dir,'CAP_states_coincidence_z.png');
end
exportgraphics(gcf,fig_path,'ContentType','image','Resolution',600)

%% now make some overview figures for centroid differences
tile_titles = {'Rest, M+RETROICOR','Rest, +PRF','Rest, +CompCor','Breathing, M+RETROICOR','Breathing, +PRF','Breathing, +CompCor'};
if z_scoring == 0
    comp_type_labels = {'vigilance','framewise displacement','respiratory flow (z-scored)','heart rate'};
elseif z_scoring == 1
    comp_type_labels = {'vigilance (z-scored)','framewise displacement (z-scored)','respiratory flow (z-scored)','heart rate (z-scored)'};
end
for m = 1:4
    centroid_diff_fig = figure('Position',[200 200 660 256]);
    t = tiledlayout(2,numel(cleanings));
    t_title = strcat('Differences associated with',{' '},comp_type_labels{m});
    title(t,t_title)
    for i = 1:2
        for c_i = 1:3
            nexttile
            
            c_limit = max(abs(centroid_diff_avgs(m,i,c_i,:,:)),[],'all');
            c_limit = 2;
            h_im = imagesc(squeeze(centroid_diff_avgs(m,i,c_i,:,:)),[-c_limit c_limit]);
            nan_spots = ~isnan(squeeze(centroid_diff_avgs(m,i,c_i,:,:)));
            set(gca, 'Color', [0.5 0.5 0.5]); % set background to gray
            set(h_im, 'AlphaData', nan_spots);
            colormap(cmap)
            c_i + 2*(i-1)
            
            ylabel('# of states')
            xticks(1:2:31)
            xlabellings = string([1:12,1:19]);
            xticklabels(xlabellings(1:2:end))
            xlabel('Brainstem          CAN-associated')
            set(gca, 'FontSize', 5)
            ax = gca;
            ax.XAxis.LabelHorizontalAlignment = 'left';
            title(tile_titles{c_i + 3*(i-1)},'FontSize',8)
            
        end
    end
    colorbar
    if z_scoring == 0
        fig_path = strcat(fig_dir,'CAP_states_assessment_centroid_diffs_',comp_types{m},'.png');
    elseif z_scoring == 1
        fig_path = strcat(fig_dir,'CAP_states_assessment_centroid_diffs_',comp_types{m},'_z.png');
    end
    exportgraphics(gcf,fig_path,'ContentType','image','Resolution',600)
end

%% make a histogram overview figure for centroid differences
centroid_diffs_all = {centroid_diffs_vig,centroid_diffs_fd,centroid_diffs_rf,centroid_diffs_hr};
figure
t = tiledlayout(numel(comp_types),2);
for m = 1:4
    for i = 1:2
        centroid_diffs_m = {};
        for c_i = 1:3
            centroid_diffs_i = [];
            centroid_diffs_i_all = centroid_diffs_all{m}{c_i+2*(i-1)};
            for ind = 1:numel(centroid_diffs_i_all)
                 centroid_diffs_i = [centroid_diffs_i;centroid_diffs_i_all{ind}];
            end            
            centroid_diffs_m{end+1} = centroid_diffs_i;
            if size(centroid_diffs_i,1) < 2
                centroid_diffs_m{end} = [centroid_diffs_i;centroid_diffs_i];
            end
        end
        nexttile
        daboxplot(centroid_diffs_m)
        hold on
        yline(0)
        title(comp_type_labels{m})
    end
end