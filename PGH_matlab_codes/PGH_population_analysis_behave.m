%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function PGH_population_analysis
function PGH_population_analysis_behave
clc; clear all; close all;
tic
subject = {'data_59d', 'data_125d'};
for counter_monkey = 1 : length(subject)
    path_data_monkey_sorted = subject{counter_monkey};
    %% params funcs
    PGH_params_funcs;
    %% Build functions
    %extract_population_data(path_data_monkey_sorted);
    %build_population_data(path_data_monkey_sorted)
    %generate_population_meta_data(path_data_monkey_sorted)
    %build_foraging_analysis(params,path_data_monkey_sorted)
    %% Plot Functions
    %plot_session_analysis(params,path_data_monkey_sorted)
    plot_foraging_analysis(params,path_data_monkey_sorted)
    toc
end
end
%% UTIL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function combine_recs
function [LICKS_ALL_DATA,SACS_ALL_DATA,BLINKS_ALL_DATA,PUPILS_ALL_DATA] = PGH_combine_recs(data_recording)
num_recording = numel(data_recording);

% LICKS_ALL_DATA
field_names_LICKS_ALL_DATA = fieldnames(data_recording(1).LICKS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_LICKS_ALL_DATA)
    field_name_LICKS_ALL_DATA = field_names_LICKS_ALL_DATA{counter_field};
    data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = horzcat(...
            data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA), ...
            data_recording(counter_recording).LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA));
    end
end
LICKS_ALL_DATA = data_cell.LICKS_ALL_DATA;

% SACS_ALL_DATA
field_names_SACS_ALL_DATA = fieldnames(data_recording(1).SACS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
    field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
    data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
            data_recording(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
    end
end
SACS_ALL_DATA = data_cell.SACS_ALL_DATA;

% BLINKS_ALL_DATA
field_names_BLINKS_ALL_DATA = fieldnames(data_recording(1).BLINKS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_BLINKS_ALL_DATA)
    field_name_BLINKS_ALL_DATA = field_names_BLINKS_ALL_DATA{counter_field};
    data_cell.BLINKS_ALL_DATA.(field_name_BLINKS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.BLINKS_ALL_DATA.(field_name_BLINKS_ALL_DATA) = horzcat(...
            data_cell.BLINKS_ALL_DATA.(field_name_BLINKS_ALL_DATA), ...
            data_recording(counter_recording).BLINKS_ALL_DATA.(field_name_BLINKS_ALL_DATA));
    end
end
BLINKS_ALL_DATA = data_cell.BLINKS_ALL_DATA;

% PUPILS_ALL_DATA
field_names_PUPILS_ALL_DATA = fieldnames(data_recording(1).PUPILS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_PUPILS_ALL_DATA)
    field_name_PUPILS_ALL_DATA = field_names_PUPILS_ALL_DATA{counter_field};
    data_cell.PUPILS_ALL_DATA.(field_name_PUPILS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.PUPILS_ALL_DATA.(field_name_PUPILS_ALL_DATA) = horzcat(...
            data_cell.PUPILS_ALL_DATA.(field_name_PUPILS_ALL_DATA), ...
            data_recording(counter_recording).PUPILS_ALL_DATA.(field_name_PUPILS_ALL_DATA));
    end
end
PUPILS_ALL_DATA = data_cell.PUPILS_ALL_DATA;
end

%% function build_variables
function [count,lick,sac] = build_variables(population_experiment_data, population_lick_data,population_sac_data,population_blink_data,population_pupil_data, path)
count.num_session =  length(population_lick_data.tag);
count.num_type = 1; % determine 1 = between harvest, 2 = within harvest trial/sac
count.num_tag_lick = 7; % 7 = consider tags 1 : 7
count.num_tag_sac = 10;

%%% compute vigor fit %%%
[FIT] = compute_vigor(population_lick_data, population_sac_data, count, path,0);
is_end_after_vigor = 0;
if is_end_after_vigor == 1
    pause;
end

for counter_session = 1 : count.num_session
    %%% build lick data %%%
    ind_str_harvest = find(population_lick_data.tag_harvest{counter_session, 1} == 1);
    ind_end_harvest = find(population_lick_data.tag_harvest{counter_session, 1} == 2);
    count.num_harvest(counter_session,1) = length(ind_str_harvest);
    count.num_trial(counter_session,1) = population_experiment_data.num_trial_sess{counter_session,1};
    count.num_lick(counter_session,1) = population_experiment_data.num_licks_sess{counter_session,1};
    count.weight(counter_session,1) = population_experiment_data.weight{counter_session,1};
    for counter_harvest = 1 : count.num_harvest(counter_session,1)
        ind_harvest_span = ind_str_harvest(counter_harvest):ind_end_harvest(counter_harvest); % span of inds from harvest str to harvest end
        lick.tag_lick{counter_session,1}{counter_harvest,:} = population_lick_data.tag{counter_session, 1}(ind_harvest_span);
        lick.tag_harvest{counter_session,1}{counter_harvest,:} = population_lick_data.tag_harvest{counter_session, 1}(ind_harvest_span);
        lick.time_onset_lick{counter_session,1}{counter_harvest,:} = population_lick_data.time_onset{counter_session, 1}(ind_harvest_span);
        lick.time_offset_lick{counter_session,1}{counter_harvest,:} = population_lick_data.time_offset{counter_session, 1}(ind_harvest_span);
        lick.tongue_dm_max{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_dm_max{counter_session, 1}(ind_harvest_span);
        lick.tongue_vm_max{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_vm_max{counter_session, 1}(ind_harvest_span);
        lick.tongue_vm_min{counter_session,1}{counter_harvest,:} = abs(population_lick_data.tongue_vm_min{counter_session, 1}(ind_harvest_span));
        lick.tongue_ang_max{counter_session,1}{counter_harvest,:} = abs(population_lick_data.tongue_ang_max{counter_session, 1}(ind_harvest_span));
        lick.ILI{counter_session,1}{counter_harvest,:} = [diff(lick.time_onset_lick{counter_session,1}{counter_harvest,:}) nan]; ...
            lick.ILI{counter_session,1}{counter_harvest,1}(1,lick.ILI{counter_session,1}{counter_harvest,1}(1,:) < 0) = nan; lick.ILI{counter_session,1}{counter_harvest,1}(1,lick.ILI{counter_session,1}{counter_harvest,1}(1,:) > 0.6) = nan;
        lick.ILR{counter_session,1}{counter_harvest,:} = 1./lick.ILI{counter_session,1}{counter_harvest,:};
        count.num_lick_harvest{counter_session,1}(counter_harvest,1) = length(lick.time_onset_lick{counter_session,1}{counter_harvest,:});
        lick.num_lick_harvest{counter_session,1}{counter_harvest,1} =  count.num_lick_harvest{counter_session,1}(counter_harvest,1); % copy for sake of figures
        lick.duration_harvest{counter_session,1}{counter_harvest,1} =  lick.time_offset_lick{counter_session,1}{counter_harvest,1}(1,end) -  lick.time_onset_lick{counter_session,1}{counter_harvest,1}(1,1);
        lick.tongue_dm{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_dm{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_vm{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_vm{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_ang{counter_session,1}{counter_harvest,1} = abs(population_lick_data.tongue_ang{counter_session,1}(:,ind_harvest_span))';
        lick.tongue_duration{counter_session,1}{counter_harvest,1} = lick.time_offset_lick{counter_session,1}{counter_harvest,:} - lick.time_onset_lick{counter_session,1}{counter_harvest,:};
        lick.tongue_tip_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_tip_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_tip_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_tip_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_mid_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_mid_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_mid_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_mid_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_r_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_r_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_r_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_r_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_l_px{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_l_px{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_l_py{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_l_py{counter_session,1}(:,ind_harvest_span)';
        lick.tongue_tip_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_tip_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_tip_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_tip_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_mid_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_mid_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_mid_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_mid_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.tongue_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.tongue_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.nose_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.nose_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.rtube_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_py_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_px_dmax{counter_session, 1}(ind_harvest_span);
        lick.ltube_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_py_dmax{counter_session, 1}(ind_harvest_span);
        rtube_r_px_onset{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_px_onset{counter_session, 1}(ind_harvest_span);
        rtube_r_py_onset{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_r_py_onset{counter_session, 1}(ind_harvest_span);
        rtube_l_px_onset{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_px_onset{counter_session, 1}(ind_harvest_span);
        rtube_l_py_onset{counter_session,1}{counter_harvest,:} = population_lick_data.rtube_l_py_onset{counter_session, 1}(ind_harvest_span);
        ltube_r_px_onset{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_px_onset{counter_session, 1}(ind_harvest_span);
        ltube_r_py_onset{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_r_py_onset{counter_session, 1}(ind_harvest_span);
        ltube_l_px_onset{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_px_onset{counter_session, 1}(ind_harvest_span);
        ltube_l_py_onset{counter_session,1}{counter_harvest,:} = population_lick_data.ltube_l_py_onset{counter_session, 1}(ind_harvest_span);

        %find if a bout is made to the right or left tube using sum of angles, then nan out the opposite reward data. if equal, completley nan out
        lick.tongue_ang_max_dir{counter_session,1}{counter_harvest,1} = population_lick_data.tongue_ang_max{counter_session,1}(:,ind_harvest_span);
        if sum(lick.tongue_ang_max_dir{counter_session,1}{counter_harvest,:} > 0) > sum(lick.tongue_ang_max_dir{counter_session,1}{counter_harvest,:} < 0)
            rew_r{counter_session,1}{counter_harvest,:} = population_lick_data.rew_capacity_r_lick_onset{counter_session, 1}(ind_harvest_span);
            rew_l{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_r_px_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_r_py_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_rtube_px{counter_session,1}{counter_harvest,:} = (rtube_r_px_onset{counter_session,1}{counter_harvest,1} + rtube_l_px_onset{counter_session,1}{counter_harvest,1})/2;
            cent_rtube_py{counter_session,1}{counter_harvest,:} = (rtube_r_py_onset{counter_session,1}{counter_harvest,1} + rtube_l_py_onset{counter_session,1}{counter_harvest,1})/2;
            cent_ltube_px{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_ltube_py{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            count.dir_harvest{counter_session,1}(counter_harvest,:) = 'r';
        elseif sum(lick.tongue_ang_max_dir{counter_session,1}{counter_harvest,:} > 0) < sum(lick.tongue_ang_max_dir{counter_session,1}{counter_harvest,:} < 0)
            rew_r{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            rew_l{counter_session,1}{counter_harvest,:} = population_lick_data.rew_capacity_l_lick_onset{counter_session, 1}(ind_harvest_span);
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_l_px_dmax{counter_session, 1}(ind_harvest_span);
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = population_lick_data.rew_l_py_dmax{counter_session, 1}(ind_harvest_span);
            cent_ltube_px{counter_session,1}{counter_harvest,:} = (ltube_r_px_onset{counter_session,1}{counter_harvest,1} + ltube_l_px_onset{counter_session,1}{counter_harvest,1})/2;
            cent_ltube_py{counter_session,1}{counter_harvest,:} = (ltube_r_py_onset{counter_session,1}{counter_harvest,1} + ltube_l_py_onset{counter_session,1}{counter_harvest,1})/2;
            cent_rtube_px{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_rtube_py{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            count.dir_harvest{counter_session,1}(counter_harvest,:) = 'l';
        else
            rew_r{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            rew_l{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_r_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_r_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_l_px_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            lick.rew_l_py_dmax{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_ltube_px{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_ltube_py{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_rtube_px{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            cent_rtube_py{counter_session,1}{counter_harvest,:} = nan(1,count.num_lick_harvest{counter_session,1}(counter_harvest,1));
            count.dir_harvest{counter_session,1}(counter_harvest,:) = 'N';
        end

        count.dist_r_tube_harvest{counter_session,1}(counter_harvest,1) = sqrt(nanmean(cent_rtube_px{counter_session,1}{counter_harvest,:})^2 + nanmean(cent_rtube_py{counter_session,1}{counter_harvest,:})^2); % as harvest series
        count.dist_l_tube_harvest{counter_session,1}(counter_harvest,1) = sqrt(nanmean(cent_ltube_px{counter_session,1}{counter_harvest,:})^2 + nanmean(cent_ltube_py{counter_session,1}{counter_harvest,:})^2); % as harvest series
        rew_r_str_harvest{counter_session,1}{counter_harvest,1} = rew_r{counter_session,1}{counter_harvest,:}(1); % as harvest series
        rew_r_end_harvest{counter_session,1}{counter_harvest,1} = rew_r{counter_session,1}{counter_harvest,:}(end); % as harvest series
        rew_l_str_harvest{counter_session,1}{counter_harvest,1} = rew_l{counter_session,1}{counter_harvest,:}(1); % as harvest series
        rew_l_end_harvest{counter_session,1}{counter_harvest,1} = rew_l{counter_session,1}{counter_harvest,:}(end); % as harvest series
        rew_r_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_r_str_harvest{counter_session,1}{counter_harvest,1} - rew_r_end_harvest{counter_session,1}{counter_harvest,1}; % as harvest series
        rew_l_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_l_str_harvest{counter_session,1}{counter_harvest,1} - rew_l_end_harvest{counter_session,1}{counter_harvest,1}; % as harvest series
        rew_r_consumed{counter_session,1}{counter_harvest,:} = [rew_r_str_harvest{counter_session,1}{counter_harvest,1} -  rew_r{counter_session,1}{counter_harvest,:}(1,2:end) nan]; % as lick series
        rew_l_consumed{counter_session,1}{counter_harvest,:} = [rew_l_str_harvest{counter_session,1}{counter_harvest,1} -  rew_l{counter_session,1}{counter_harvest,:}(1,2:end) nan]; % as lick series

        % combine food tube data
        if ~isnan(rew_r_str_harvest{counter_session,1}{counter_harvest,1})
            lick.rew{counter_session,1}{counter_harvest,:} = rew_r{counter_session,1}{counter_harvest,:};
            lick.rew_str_harvest{counter_session,1}{counter_harvest,1} = rew_r_str_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_end_harvest{counter_session,1}{counter_harvest,1} = rew_r_end_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_r_consumed_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_r_consumed{counter_session,1}{counter_harvest,:};
            count.dist_tube_harvest{counter_session,1}(counter_harvest,1) =  count.dist_r_tube_harvest{counter_session,1}(counter_harvest,1);
        elseif ~isnan(rew_l_str_harvest{counter_session,1}{counter_harvest,1})
            lick.rew{counter_session,1}{counter_harvest,:} = rew_l{counter_session,1}{counter_harvest,:};
            lick.rew_str_harvest{counter_session,1}{counter_harvest,1} = rew_l_str_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_end_harvest{counter_session,1}{counter_harvest,1} = rew_l_end_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_l_consumed_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_l_consumed{counter_session,1}{counter_harvest,:};
            count.dist_tube_harvest{counter_session,1}(counter_harvest,1) =  count.dist_l_tube_harvest{counter_session,1}(counter_harvest,1);
        else % use r since it will be nan anyway
            lick.rew{counter_session,1}{counter_harvest,:} = rew_r{counter_session,1}{counter_harvest,:};
            lick.rew_str_harvest{counter_session,1}{counter_harvest,1} = rew_r_str_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_end_harvest{counter_session,1}{counter_harvest,1} = rew_r_end_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = rew_r_consumed_harvest{counter_session,1}{counter_harvest,1};
            lick.rew_consumed{counter_session,1}{counter_harvest,:} = rew_r_consumed{counter_session,1}{counter_harvest,:};
            count.dist_tube_harvest{counter_session,1}(counter_harvest,1) =  count.dist_r_tube_harvest{counter_session,1}(counter_harvest,1);
        end

        lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1}(lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1}<0) = nan;
        lick.rew_consumed{counter_session,1}{counter_harvest,1}(lick.rew_consumed{counter_session,1}{counter_harvest,1}<0) = nan;

        lick.rew_consumed_rate_harvest{counter_session,1}{counter_harvest,1} = lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1}/count.num_lick_harvest{counter_session,1}(counter_harvest,1);
        lick.rew_consumed_rate_harvest{counter_session,1}{counter_harvest,1}(lick.rew_consumed_rate_harvest{counter_session,1}{counter_harvest,1} > 0.1) = nan;
        lick.rew_consumed_rate{counter_session,1}{counter_harvest,1} = [nan diff(lick.rew_consumed{counter_session,1}{counter_harvest,1})];
        lick.rew_consumed_rate{counter_session,1}{counter_harvest,1}(lick.rew_consumed_rate{counter_session,1}{counter_harvest,1} > 0.6 | lick.rew_consumed_rate{counter_session,1}{counter_harvest,1} < -0.6 ) = nan;

        %%% load and compute vigor of pro and ret %%%
        load([path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_FIT.mat' ]);
        for counter_lick = 1 : count.num_lick_harvest{counter_session,1}(counter_harvest,1)
            if (lick.tag_lick{counter_session,1}{counter_harvest,:}(1,counter_lick) == 1)
                fit_pro = FIT.fit_lick_groom_pro ;
                fit_ret = FIT.fit_lick_groom_ret ;
            else
                fit_pro = FIT.fit_lick_rew_pro ;
                fit_ret = FIT.fit_lick_rew_ret ;
            end
            lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,counter_lick) = lick.tongue_vm_max{counter_session,1}{counter_harvest,:}(1,counter_lick) ...
                /fit_pro(lick.tongue_dm_max{counter_session,1}{counter_harvest,:}(1,counter_lick));
            lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,counter_lick) = lick.tongue_vm_min{counter_session,1}{counter_harvest,:}(1,counter_lick) ...
                /fit_ret(lick.tongue_dm_max{counter_session,1}{counter_harvest,:}(1,counter_lick));
        end

        lick.tongue_vigor_pro{counter_session,1}{counter_harvest,1}(1, lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.tongue_vigor_pro{counter_session,1}{counter_harvest,:}(1,:) < 0.5) = nan;
        lick.tongue_vigor_ret{counter_session,1}{counter_harvest,1}(1, lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,:) > 1.5 | lick.tongue_vigor_ret{counter_session,1}{counter_harvest,:}(1,:) < 0.5) = nan;

        %%% compute reward gained during work %%%
        if counter_harvest == 1
            sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan ; % as harvest series
        else
            if ~isnan(lick.rew_str_harvest{counter_session,1}{counter_harvest,1}) && ~isnan(lick.rew_end_harvest{counter_session,1}{counter_harvest-1,1})
                sac.rew_gained_work{counter_session,1}{counter_harvest,1} = lick.rew_str_harvest{counter_session,1}{counter_harvest,1} - lick.rew_end_harvest{counter_session,1}{counter_harvest-1,1};
                if  sac.rew_gained_work{counter_session,1}{counter_harvest,1} < 0 || sac.rew_gained_work{counter_session,1}{counter_harvest,1} > 1.25
                    sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan;
                end
            else
                sac.rew_gained_work{counter_session,1}{counter_harvest,1} = nan;
            end
        end

        %%% filter rew %%%
        lick.rew_consumed{counter_session,1}{counter_harvest,1}(1, lick.rew_consumed{counter_session,1}{counter_harvest,:}(1,:) > 1.25 | lick.rew_consumed{counter_session,1}{counter_harvest,:}(1,:) < 0) = nan;
        lick.rew{counter_session,1}{counter_harvest,1}(1, lick.rew{counter_session,1}{counter_harvest,:}(1,:) > 1.25 | lick.rew{counter_session,1}{counter_harvest,:}(1,:) < 0) = nan;

        if (lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} > 1.25 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < 0)
            lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} = nan;
        end
        if (lick.rew_str_harvest{counter_session,1}{counter_harvest,1} > 1.25 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < 0 )
            lick.rew_str_harvest{counter_session,1}{counter_harvest,1} = nan;
        end
        if (lick.rew_end_harvest{counter_session,1}{counter_harvest,1} > 1.25 || lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1} < 0 )
            lick.rew_end_harvest{counter_session,1}{counter_harvest,1} = nan;
        end
    end

    % take average of dist_tube r and l within a session and add to count
    count.dist_r_tube_harvest_mean(counter_session,1) =  nanmean(count.dist_r_tube_harvest{counter_session,1});
    count.dist_l_tube_harvest_mean(counter_session,1) =  nanmean(count.dist_l_tube_harvest{counter_session,1});

    % combine r and l using means of dist_tube r and l
    count.dist_tube_harvest_mean{counter_session,1} =  count.dist_r_tube_harvest{counter_session,1};
    count.dist_tube_harvest_mean{counter_session,1}(~isnan(count.dist_tube_harvest_mean{counter_session,1})) = count.dist_r_tube_harvest_mean(counter_session,1);
    count.dist_tube_harvest_mean{counter_session,1}(isnan(count.dist_tube_harvest_mean{counter_session,1})) = count.dist_l_tube_harvest_mean(counter_session,1);

    %%% build eye data %%%
    % shift eye data to align eyelink data with video data using sample_diff
    sample_diff = population_experiment_data.sample_diff{counter_session, 1};
    num_sac_rec{counter_session, 1} = population_experiment_data.num_sacs_sess_rec{counter_session, 1};
    exp_start_time   = population_experiment_data.exp_start_time{counter_session, 1};
    num_trial_rec{counter_session, 1} = population_experiment_data.num_trial_sess_rec{counter_session, 1};
    num_harvest_rec{counter_session,1} = population_experiment_data.num_harvests_sess_rec{counter_session, 1};
    num_blink_rec{counter_session,1} = population_experiment_data.num_blink_sess_rec{counter_session, 1};
    length_pupil_rec{counter_session,1} = population_experiment_data.length_pupil_sess_rec{counter_session, 1};
    time_start_trial_{counter_session,1} = population_sac_data.time_start{counter_session,1};
    time_end_trial_{counter_session,1} = population_sac_data.time_end{counter_session,1};
    tag_sac_{counter_session,1} = population_sac_data.tag{counter_session,1};
    time_onset_sac_{counter_session,1} = population_sac_data.time_onset{counter_session,1};
    time_onset_fix_{counter_session,1} = population_sac_data.time_fix_onset{counter_session,1};
    time_offset_fix_{counter_session,1} = population_sac_data.time_fix_offset{counter_session,1};
    eye_vm_{counter_session,1} = population_sac_data.eye_r_vm{counter_session,1};
    eye_ang_{counter_session,1} = population_sac_data.eye_r_ang{counter_session,1};
    eye_vm_max_{counter_session,1} = max(population_sac_data.eye_r_vm{counter_session,1});
    eye_dm_max_{counter_session,1} = population_sac_data.eye_r_amp_m{counter_session,1};
    trial_num_{counter_session,1} = population_sac_data.trial_num{counter_session,1};
    num_trial_attempt_rec_{counter_session, 1} = population_sac_data.num_trial_attempt{counter_session, 1};
    reaction_{counter_session,1} = population_sac_data.reaction{counter_session,1};
    eye_px_onset_{counter_session,1} = population_sac_data.eye_r_px_onset{counter_session,1};
    eye_py_onset_{counter_session,1} = population_sac_data.eye_r_py_onset{counter_session,1};
    tgt_px_onset_{counter_session,1} = population_sac_data.tgt_px_onset{counter_session,1};
    tgt_py_onset_{counter_session,1} = population_sac_data.tgt_py_onset{counter_session,1};
    eye_px_offset_{counter_session,1} = population_sac_data.eye_r_px_offset{counter_session,1};
    eye_py_offset_{counter_session,1} = population_sac_data.eye_r_py_offset{counter_session,1};
    tgt_px_offset_{counter_session,1} = population_sac_data.tgt_px_offset{counter_session,1};
    tgt_py_offset_{counter_session,1} = population_sac_data.tgt_py_offset{counter_session,1};
    eye_px_before_{counter_session,1} = population_sac_data.x_fix_before{counter_session,1};
    eye_py_before_{counter_session,1} = population_sac_data.y_fix_before{counter_session,1};
    eye_px_after_{counter_session,1} = population_sac_data.x_fix_after{counter_session,1};
    eye_py_after_{counter_session,1} = population_sac_data.y_fix_after{counter_session,1};
    validity_sac_{counter_session,1} = population_sac_data.validity{counter_session,1};
    validity_fix_{counter_session,1} = population_sac_data.fix_validity{counter_session,1};
    time_onset_blink_{counter_session,1} = population_blink_data.blink_onset{counter_session, 1};
    time_offset_blink_{counter_session,1} = population_blink_data.blink_offset{counter_session, 1};
    time_1K_{counter_session,1} = population_pupil_data.time_1K{counter_session,1};
    validity_pupil_{counter_session,1} = population_pupil_data.validity{counter_session,1};
    pupil_area_raw_{counter_session,1} = population_pupil_data.pupil_area{counter_session,1};
    % filter pupil area
    pupil_area_{counter_session,1} = filter_pupil(pupil_area_raw_{counter_session,1},counter_session,path);
    % z score pupil area
    pupil_area_{counter_session,1} = (pupil_area_{counter_session,1}-nanmean(pupil_area_{counter_session,1}))/nanstd(pupil_area_{counter_session,1});

    shift_sac = 1;
    shift_trial = 1;
    shift_blink = 1;
    shift_pupil = 1;

    for counter_rec = 1 :  length(num_trial_rec{counter_session})
        time_start_trial__{counter_session,1}{1,counter_rec} = (time_start_trial_{counter_session,1}(1,shift_trial : num_trial_rec{counter_session, 1}(counter_rec) + shift_trial - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_end_trial__{counter_session,1}{1,counter_rec} = (time_end_trial_{counter_session,1}(1,shift_trial : num_trial_rec{counter_session, 1}(counter_rec) + shift_trial - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_onset_sac__{counter_session,1}{1,counter_rec} = (time_onset_sac_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_onset_fix__{counter_session,1}{1,counter_rec} = (time_onset_fix_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_offset_fix__{counter_session,1}{1,counter_rec} = (time_offset_fix_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        tag_sac__{counter_session,1}{1,counter_rec} = (tag_sac_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_vm__{counter_session,1}{1,counter_rec} = (eye_vm_{counter_session,1}(:,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_ang__{counter_session,1}{1,counter_rec} = (eye_ang_{counter_session,1}(:,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_vm_max__{counter_session,1}{1,counter_rec} = (eye_vm_max_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_dm_max__{counter_session,1}{1,counter_rec} = (eye_dm_max_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        trial_num__{counter_session,1}{1,counter_rec} = (trial_num_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        num_trial_attempt__{counter_session,1}{1,counter_rec} = (num_trial_attempt_rec_{counter_session,1}(1,shift_trial : num_trial_rec{counter_session, 1}(counter_rec) + shift_trial - 1 ));
        reaction__{counter_session,1}{1,counter_rec} = (reaction_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 )); ...
            reaction__{counter_session,1}{1,counter_rec}(reaction__{counter_session,1}{1,counter_rec} <0 | reaction__{counter_session,1}{1,counter_rec} > 600) = nan;
        eye_px_onset__{counter_session,1}{1,counter_rec} = (eye_px_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_onset__{counter_session,1}{1,counter_rec} = (eye_py_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_px_onset__{counter_session,1}{1,counter_rec} = (tgt_px_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_py_onset__{counter_session,1}{1,counter_rec} = (tgt_py_onset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_px_offset__{counter_session,1}{1,counter_rec} = (eye_px_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_offset__{counter_session,1}{1,counter_rec} = (eye_py_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_px_offset__{counter_session,1}{1,counter_rec} = (tgt_px_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        tgt_py_offset__{counter_session,1}{1,counter_rec} = (tgt_py_offset_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_px_before__{counter_session,1}{1,counter_rec} = (eye_px_before_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_before__{counter_session,1}{1,counter_rec} = (eye_py_before_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_px_after__{counter_session,1}{1,counter_rec} = (eye_px_after_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        eye_py_after__{counter_session,1}{1,counter_rec} = (eye_py_after_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        validity_sac__{counter_session,1}{1,counter_rec} = (validity_sac_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));
        validity_fix__{counter_session,1}{1,counter_rec} = (validity_fix_{counter_session,1}(1,shift_sac : num_sac_rec{counter_session, 1}(counter_rec) + shift_sac - 1 ));

        %         time_onset_blink__{counter_session,1}{1,counter_rec} = (time_onset_blink_{counter_session,1}(1,shift_blink : num_blink_rec{counter_session, 1}(counter_rec) + shift_blink - 1 ) ...
        %             - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        %         time_offset_blink__{counter_session,1}{1,counter_rec} = (time_offset_blink_{counter_session,1}(1,shift_blink : num_blink_rec{counter_session, 1}(counter_rec) + shift_blink - 1 ) ...
        %             - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        time_1K__{counter_session,1}{1,counter_rec} = (time_1K_{counter_session,1}(1,shift_pupil : length_pupil_rec{counter_session, 1}(counter_rec) + shift_pupil - 1 ) ...
            - exp_start_time(counter_rec) + sample_diff(counter_rec)/1000);
        pupil_area__{counter_session,1}{1,counter_rec} = (pupil_area_{counter_session,1}(1,shift_pupil : length_pupil_rec{counter_session, 1}(counter_rec) + shift_pupil - 1 ));
        validity_pupil__{counter_session,1}{1,counter_rec} = (validity_pupil_{counter_session,1}(1,shift_pupil : length_pupil_rec{counter_session, 1}(counter_rec) + shift_pupil - 1 ));

        shift_sac = shift_sac + num_sac_rec{counter_session, 1}(counter_rec);
        shift_trial = shift_trial + num_trial_rec{counter_session, 1}(counter_rec);
        shift_blink = shift_blink + num_blink_rec{counter_session, 1}(counter_rec);
        shift_pupil = shift_pupil + length_pupil_rec{counter_session, 1}(counter_rec);

        for counter_harvest_rec = 1 : num_harvest_rec{counter_session,1}(counter_rec)
            shift = 0;
            if counter_rec > 1
                counter_harvest_rec = counter_harvest_rec + sum(num_harvest_rec{counter_session,1}(1 : counter_rec - 1));
                shift = sum(num_harvest_rec{counter_session,1}(1 : counter_rec - 1));
            end

            % detect trials and saccs occuring either pre or within bout
            if counter_harvest_rec == 1
                is_pre_bout_trial = time_end_trial__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);
                is_pre_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);
                %is_pre_bout_blink = time_onset_blink__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);
                is_pre_bout_pupil = time_1K__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1);

            elseif counter_harvest_rec > 1
                is_pre_bout_trial = time_end_trial__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                    & time_end_trial__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
                is_pre_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                    & time_onset_sac__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
                % is_pre_bout_blink = time_onset_blink__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                % & time_onset_blink__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
                is_pre_bout_pupil = time_1K__{counter_session,1}{counter_rec} < lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                    & time_1K__{counter_session,1}{counter_rec} > lick.time_onset_lick{counter_session,1}{counter_harvest_rec-1,:}(1,end);
            end
            is_in_bout_trial = time_end_trial__{counter_session,1}{counter_rec} >= lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                & time_end_trial__{counter_session,1}{counter_rec} <= lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);
            is_in_bout_sac = time_onset_sac__{counter_session,1}{counter_rec} >= lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                & time_onset_sac__{counter_session,1}{counter_rec} <= lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);
            % is_in_bout_blink = time_onset_blink__{counter_session,1}{counter_rec} >= lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
            % & time_onset_blink__{counter_session,1}{counter_rec} <= lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);
            is_in_bout_pupil = time_1K__{counter_session,1}{counter_rec} >= lick.time_onset_lick{counter_session,1}{counter_harvest_rec,:}(1,1) ...
                & time_1K__{counter_session,1}{counter_rec} <= lick.time_offset_lick{counter_session,1}{counter_harvest_rec,:}(1,end);

            sac.time_start_trial{counter_session,1}{counter_harvest_rec,:} = time_start_trial__{counter_session,1}{counter_rec}(1,is_pre_bout_trial);
            sac.time_end_trial{counter_session,1}{counter_harvest_rec,:} = time_end_trial__{counter_session,1}{counter_rec}(1,is_pre_bout_trial);

            count.num_trial_work{counter_session,1}{counter_harvest_rec,1} = length(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:});
            sac.num_trial_work{counter_session,1}{counter_harvest_rec,1} = count.num_trial_work{counter_session,1}{counter_harvest_rec,1};
            sac.num_trial_work{counter_session,1}{counter_harvest_rec,1}(sac.num_trial_work{counter_session,1}{counter_harvest_rec,1}>15) = 5;
            count.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1} = num_trial_attempt__{counter_session,1}{counter_rec}(1,is_pre_bout_trial);
            sac.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1} = sum(count.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1});
            sac.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1}(sac.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1}>40) = 8;

            % handle situations where no trials occur
            if isempty(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:})
                sac.ITI{counter_session,1}{counter_harvest_rec,:} = [];
                sac.ITR{counter_session,1}{counter_harvest_rec,:} = [];
                count.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1} = 0;
                sac.num_trial_attempt_work{counter_session,1}{counter_harvest_rec,1} = 0;
            else
                sac.ITI{counter_session,1}{counter_harvest_rec,:} = [diff(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:}) nan]; sac.ITI{counter_session,1}{counter_harvest_rec,:}(sac.ITI{counter_session,1}{counter_harvest_rec,:} <= 0 | sac.ITI{counter_session,1}{counter_harvest_rec,:} > 60) = nan;
                sac.ITR{counter_session,1}{counter_harvest_rec,:} = 1./sac.ITI{counter_session,1}{counter_harvest_rec,:};
            end

            % specify saccade type to analyze
            tag_146 = [1 4 6];
            tag_10 = [10];

            is_tag_146 = ismember(tag_sac__{counter_session,1}{counter_rec}, tag_146);
            if tag_sac__{counter_session,1}{counter_rec}(find(is_pre_bout_sac == 1, 1, 'first')) == 4
                is_tag_146(find(is_pre_bout_sac == 1, 1, 'first')) = 0;
            end

            is_tag_10 = ismember(tag_sac__{counter_session,1}{counter_rec}, tag_10);

            % specify validity
            is_valid_sac = validity_sac__{counter_session,1}{counter_rec} == 1;
            is_valid_fix = validity_fix__{counter_session,1}{counter_rec} == 1;
            is_valid_pupil = validity_pupil__{counter_session,1}{counter_rec} == 1; % add pupil area filter

            sac.tag_sac{counter_session,1}{counter_harvest_rec,:} = tag_sac__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.trial_num{counter_session,1}{counter_harvest_rec,:} = trial_num__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:} = time_onset_sac__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.time_onset_fix_h{counter_session,1}{counter_harvest_rec,:} = time_onset_fix__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.time_offset_fix_h{counter_session,1}{counter_harvest_rec,:} = time_offset_fix__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.eye_duration_fix_h{counter_session,1}{counter_harvest_rec,:} =  sac.time_offset_fix_h{counter_session,1}{counter_harvest_rec,:} -  sac.time_onset_fix_h{counter_session,1}{counter_harvest_rec,:};
            sac.eye_vm_max_w{counter_session,1}{counter_harvest_rec,:} = eye_vm_max__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_vm_max_h{counter_session,1}{counter_harvest_rec,:} = eye_vm_max__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag_10 & is_valid_sac);
            sac.eye_ang{counter_session,1}{counter_harvest_rec,:} = eye_ang__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_dm_max_w{counter_session,1}{counter_harvest_rec,:} = eye_dm_max__{counter_session,1}{counter_rec}(1,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_dm_max_h{counter_session,1}{counter_harvest_rec,:} = eye_dm_max__{counter_session,1}{counter_rec}(1,is_in_bout_sac & is_tag_10 & is_valid_sac);
            sac.eye_vm{counter_session,1}{counter_harvest_rec,:} = eye_vm__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.reaction{counter_session,1}{counter_harvest_rec,:} = reaction__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_px_onset{counter_session,1}{counter_harvest_rec,:} = eye_px_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_py_onset{counter_session,1}{counter_harvest_rec,:} = eye_py_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.tgt_px_onset{counter_session,1}{counter_harvest_rec,:} = tgt_px_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.tgt_py_onset{counter_session,1}{counter_harvest_rec,:} = tgt_py_onset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_px_offset{counter_session,1}{counter_harvest_rec,:} = eye_px_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_py_offset{counter_session,1}{counter_harvest_rec,:} = eye_py_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.tgt_px_offset{counter_session,1}{counter_harvest_rec,:} = tgt_px_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.tgt_py_offset{counter_session,1}{counter_harvest_rec,:} = tgt_py_offset__{counter_session,1}{counter_rec}(:,is_pre_bout_sac & is_tag_146 & is_valid_sac);
            sac.eye_px_before_h{counter_session,1}{counter_harvest_rec,:} = eye_px_before__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.eye_py_before_h{counter_session,1}{counter_harvest_rec,:} = eye_py_before__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.eye_px_after_h{counter_session,1}{counter_harvest_rec,:} = eye_px_after__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.eye_py_after_h{counter_session,1}{counter_harvest_rec,:} = eye_py_after__{counter_session,1}{counter_rec}(:,is_in_bout_sac & is_tag_10 & is_valid_fix);
            sac.eye_vigor_w{counter_session,1}{counter_harvest_rec,:} = sac.eye_vm_max_w{counter_session,1}{counter_harvest_rec,:} ./ ...
                FIT.fit_sac_rel(sac.eye_dm_max_w{counter_session,1}{counter_harvest_rec,:})';
            sac.eye_vigor_h{counter_session,1}{counter_harvest_rec,:} = sac.eye_vm_max_h{counter_session,1}{counter_harvest_rec,:} ./ ...
                FIT.fit_sac_irr(sac.eye_dm_max_h{counter_session,1}{counter_harvest_rec,:})';
            sac.eye_vigor_w{counter_session,1}{counter_harvest_rec,:}(1, sac.eye_vigor_w{counter_session,1}{counter_harvest_rec,:}(1,:) > 1.5 | sac.eye_vigor_w{counter_session,1}{counter_harvest_rec,:}(1,:) < 0.5) = nan;
            sac.eye_vigor_h{counter_session,1}{counter_harvest_rec,:}(1, sac.eye_vigor_h{counter_session,1}{counter_harvest_rec,:}(1,:) > 1.5 | sac.eye_vigor_h{counter_session,1}{counter_harvest_rec,:}(1,:) < 0.5) = nan;
            count.num_sac_work{counter_session,1}(counter_harvest_rec,1) = length(sac.eye_vigor_w{counter_session,1}{counter_harvest_rec,:});
            count.num_sac_harvest{counter_session,1}(counter_harvest_rec,1) = length(sac.eye_vigor_h{counter_session,1}{counter_harvest_rec,:});
            count.num_fix_harvest{counter_session,1}(counter_harvest_rec,1) = length(sac.time_onset_fix_h{counter_session,1}{counter_harvest_rec,:});

            %sac.time_onset_blink{counter_session,1}{counter_harvest_rec,:} = time_onset_blink__{counter_session,1}{counter_rec}(:,is_pre_bout_blink);
            %sac.time_offset_blink{counter_session,1}{counter_harvest_rec,:} = time_offset_blink__{counter_session,1}{counter_rec}(:,is_pre_bout_blink);
            %sac.blink_duration{counter_session,1}{counter_harvest_rec,:} = sac.time_offset_blink{counter_session,1}{counter_harvest_rec,:} - sac.time_onset_blink{counter_session,1}{counter_harvest_rec,:};

            time_1K_work_{counter_session,1}{counter_harvest_rec,:} = time_1K__{counter_session,1}{counter_rec}(:,is_pre_bout_pupil & is_valid_pupil);
            time_1K_harvest_{counter_session,1}{counter_harvest_rec,:} = time_1K__{counter_session,1}{counter_rec}(:,is_in_bout_pupil & is_valid_pupil);
            pupil_area_work_{counter_session,1}{counter_harvest_rec,:} = pupil_area__{counter_session,1}{counter_rec}(:,is_pre_bout_pupil & is_valid_pupil);
            pupil_area_harvest_{counter_session,1}{counter_harvest_rec,:} = pupil_area__{counter_session,1}{counter_rec}(:,is_in_bout_pupil & is_valid_pupil);

            if isempty(sac.time_end_trial{counter_session,1}{counter_harvest_rec,:}) || isempty(sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:})
                sac.duration_work{counter_session,1}{counter_harvest_rec,:} = 0;
            else
                tags = [1 2 6 7 8];
                is_tag = ismember(sac.tag_sac{counter_session,1}{counter_harvest_rec,:},tags);
                ind_tag = find(is_tag == 1, 1, 'first');
                if isempty(ind_tag)
                    sac.duration_work{counter_session,1}{counter_harvest_rec,:} = 0;
                else
                    sac.duration_work{counter_session,1}{counter_harvest_rec,:} = sac.time_end_trial{counter_session,1}{counter_harvest_rec,:}(end) - sac.time_onset_sac{counter_session,1}{counter_harvest_rec,:}(ind_tag);
                end
                sac.duration_work{counter_session,1}{counter_harvest_rec,:}(sac.duration_work{counter_session,1}{counter_harvest_rec,:} < 0,1) = 0;
                %                 sac.duration_work{counter_session,1}{counter_harvest_rec,:}(sac.duration_work{counter_session,1}{counter_harvest_rec,:} > 60,1) = nan;
            end

            for counter_type = 1 : count.num_type
                % handle situations where no saccades occur
                if isempty(sac.time_onset_sac{counter_session,counter_type}{counter_harvest_rec,:})
                    sac.ISI{counter_session,counter_type}{counter_harvest_rec,:} = [];
                    sac.ISR{counter_session,counter_type}{counter_harvest_rec,:} = [];
                else
                    sac.ISI{counter_session,counter_type}{counter_harvest_rec,:} = [diff(sac.time_onset_sac{counter_session,counter_type}{counter_harvest_rec,:}) nan];
                    sac.ISR{counter_session,counter_type}{counter_harvest_rec,:} = 1./sac.ISI{counter_session,counter_type}{counter_harvest_rec,:};
                end
            end
        end
    end

    % calculate reward acquired rate
    for counter_harvest = 1 : count.num_harvest(counter_session)
        lick.rew_acquired_rate_harvest{counter_session,1}{counter_harvest,1} = lick.rew_consumed_harvest{counter_session,1}{counter_harvest,1}/(lick.duration_harvest{counter_session,1}{counter_harvest,1} + sac.duration_work{counter_session,1}{counter_harvest,1});
        lick.rew_acquired_rate_harvest{counter_session,1}{counter_harvest,1}(lick.rew_acquired_rate_harvest{counter_session,1}{counter_harvest,1} > 0.1) = nan;
    end

    % can be used for padding
    %     max_num_lick_harvest_sess(counter_session,1) = max(cell2mat(count.num_lick_harvest{counter_session,1}));
    %     max_num_trial_work_sess(counter_session,1) = max(cell2mat(count.num_trial_work{counter_session,1}));
    %     max_num_trial_work_sess(counter_session,2) = max(cell2mat(count.num_trial_work{counter_session,2}));
    %     max_num_sac_work_sess(counter_session,1) = max(cell2mat(count.num_sac_work{counter_session,1}));
    %     max_num_sac_work_sess(counter_session,2) = max(cell2mat(count.num_sac_work{counter_session,2}));

end

%%% build pupil data %%%
[sac.pupil_area,sac.pupil_velocity, lick.pupil_area, lick.pupil_velocity] = build_pupil(pupil_area_work_, pupil_area_harvest_,time_1K_work_,time_1K_harvest_,lick.time_onset_lick,sac.time_onset_sac,count);

%%% set max number for later padding %%%
count.max_num_harvest = max(count.num_harvest);
count.max_num_lick_harvest= 30;
count.max_num_trial_work(1,1) = 10;
count.max_num_sac_work(1,1) = 10;

end
%% function build_model
function [lick] = build_model(lick, sac,count)
for counter_session = 1 : count.num_session
    for counter_harvest = 1 : count.num_harvest(counter_session,1)
        for counter_lick = 1 : count.num_lick_harvest{counter_session,1}(counter_harvest,1)
            % define variables
            a = 60 - count.weight(counter_session, 1)/10;
            bs = sac.num_trial_work{counter_session, 1}{counter_harvest, 1}/sac.num_trial_attempt_work{counter_session, 1}{counter_harvest, 1};
            bl = sum(ismember(lick.tag_lick{counter_session,1}{counter_harvest,1}, [2 4 6]))/sum(ismember(lick.tag_lick{counter_session,1}{counter_harvest,1}, [2:7]));
            ns = sac.num_trial_work{counter_session, 1}{counter_harvest, 1};
            nl = counter_lick;
            ts = nanmean(sac.ITI{counter_session, 1}{counter_harvest, 1});
            tl = (lick.ILI{counter_session, 1}{counter_harvest, 1}(1,counter_lick));
            d = count.dist_tube_harvest{counter_session, 1}(counter_harvest,1)/10;
            k = 1;
            cs = 0.5;
            cl = d^2/tl + k*tl;

            % model
            lick.J{counter_session, 1}{counter_harvest, 1}(1,counter_lick) = (a*bs*ns*(1-1/(1+bl*nl))-nl*cl-ns^2*cs)/(nl*tl + ns*ts);

            % metabolic cost of licking
            lick.cl{counter_session, 1}{counter_harvest, 1}(1,counter_lick) = cl;
        end
    end
end
end
%% function endpoint_error_mag
function [data] = endpoint_error_mag(data, count)
fieldname = fieldnames(data);
for counter_session = 1 : count.num_session
    for counter_type = 1 : count.num_type
        if counter_type == 1
            for counter_fieldname = 1 : length(fieldname)
                % check if lick or sac related data
                if sum(contains(string(fieldname),'_lick')) > 0
                    for counter_harvest = 1 : count.num_harvest(counter_session)
                        for counter_lick = 1 : count.num_lick_harvest{counter_session,1}(counter_harvest,1)

                            % tip
                            err_tongue_tip_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_tip_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_tip_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_tip_rew_ = sqrt(nansum([err_tongue_tip_rew_r_px_; err_tongue_tip_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_tip_rew_r_py_; err_tongue_tip_rew_l_py_]).^2);
                            % mid
                            err_tongue_mid_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_mid_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_mid_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_mid_rew_ = sqrt(nansum([err_tongue_mid_rew_r_px_; err_tongue_mid_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_mid_rew_r_py_; err_tongue_mid_rew_l_py_]).^2);
                            % r
                            err_tongue_r_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_r_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_r_rew_ = sqrt(nansum([err_tongue_r_rew_r_px_; err_tongue_r_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_r_rew_r_py_; err_tongue_r_rew_l_py_]).^2);
                            % l
                            err_tongue_l_rew_r_px_ = (data.rew_r_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_r_py_ = (data.rew_r_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_l_px_ = ( data.rew_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_px_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            err_tongue_l_rew_l_py_ = (data.rew_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                - data.tongue_l_py_dmax{counter_session,1}{counter_harvest,1}(1,counter_lick)) ;
                            dist_err_tongue_l_rew_ = sqrt(nansum([err_tongue_l_rew_r_px_; err_tongue_l_rew_l_px_]).^2 ...
                                +  nansum([err_tongue_l_rew_r_py_; err_tongue_l_rew_l_py_]).^2);

                            dist_err_tongue_rew_ = [dist_err_tongue_tip_rew_; dist_err_tongue_mid_rew_; dist_err_tongue_r_rew_; dist_err_tongue_l_rew_];
                            data.dist_err_tongue_rew{counter_session,1}{counter_harvest,1}(1,counter_lick) = min(dist_err_tongue_rew_);

                            if (min(dist_err_tongue_rew_) == dist_err_tongue_tip_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_mid_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_mid_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_r_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_r_rew_l_py_;
                            elseif (min(dist_err_tongue_rew_) == dist_err_tongue_l_rew_)
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_l_rew_l_py_;
                            else
                                err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_px_;
                                err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_r_py_;
                                err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_px_;
                                err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_tongue_tip_rew_l_py_;
                            end

                            theta{counter_session,1}{counter_harvest,1}(1,counter_lick) = data.tongue_ang_max_dir{counter_session,1}{counter_harvest,1}(1,counter_lick);

                            R = [cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) -sind(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) ; ...
                                sind(theta{counter_session,1}{counter_harvest,1}(1,counter_lick)) cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_lick))];

                            err_rot = R * [err_tongue_rew_r_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) err_tongue_rew_l_px_{counter_session,1}{counter_harvest,1}(1,counter_lick) ...
                                ; err_tongue_rew_r_py_{counter_session,1}{counter_harvest,1}(1,counter_lick) err_tongue_rew_l_py_{counter_session,1}{counter_harvest,1}(1,counter_lick)];

                            if(~isnan(err_rot(1)))
                                data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(1);
                                data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(2);
                            elseif(~isnan(err_rot(3)))
                                data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(3);
                                data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = err_rot(4);
                            else
                                data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = nan;
                                data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick) = nan;
                            end

                            %data.dist_err_tongue_rew{counter_session,1}{counter_harvest,1}(1,counter_lick) = sqrt(data.err_tongue_rew_px_rot{counter_session,1}{counter_harvest,1}(1,counter_lick).^2 +  data.err_tongue_rew_py_rot{counter_session,1}{counter_harvest,1}(1,counter_lick).^2);
                            %  data. det_var_err_tongue_rew_{counter_session,1}{counter_harvest,1}(1,counter_lick) = det(nancov(data.err_tongue_rew_px_rot_{counter_session,1}{counter_harvest,1}(1,counter_lick), data.err_tongue_rew_py_rot_{counter_session,1}{counter_harvest,1}(1,counter_lick)));
                        end
                    end
                elseif sum(contains(string(fieldname),'_sac')) > 0
                    for counter_harvest = 1 : count.num_harvest(counter_session)
                        % initialize to handle situations where harvest has
                        % no saccades
                        %                         data.err_eye_tgt_px_rot{counter_session,1}(counter_harvest,1 = nan(size(count.num_sac_work{counter_harvest,1},1),1);
                        %                         data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1} =nan(1,count.num_sac_work(counter_harvest));
                        %                         data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1} = nan(1,count.num_sac_work(counter_harvest));
                        for counter_sac = 1 : count.num_sac_work{counter_session,1}(counter_harvest,1)
                            % compute x and y of error vector after rotation
                            err_eye_tgt_px{counter_session,1}{counter_harvest,1}(1,counter_sac) = (data.tgt_px_offset{counter_session,1}{counter_harvest,1}(1,counter_sac) ...
                                - data.eye_px_offset{counter_session,1}{counter_harvest,1}(1,counter_sac)) ;
                            err_eye_tgt_py{counter_session,1}{counter_harvest,1}(1,counter_sac) = (data.tgt_py_offset{counter_session,1}{counter_harvest,1}(1,counter_sac) ...
                                - data.eye_py_offset{counter_session,1}{counter_harvest,1}(1,counter_sac)) ;

                            theta{counter_session,1}{counter_harvest,1}(1,counter_sac) = data.eye_ang{counter_session,1}{counter_harvest,1}(1,counter_sac);

                            R = [cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) -sind(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) ; ...
                                sind(theta{counter_session,1}{counter_harvest,1}(1,counter_sac)) cosd(theta{counter_session,1}{counter_harvest,1}(1,counter_sac))];

                            err_rot = R * [err_eye_tgt_px{counter_session,1}{counter_harvest,1}(1,counter_sac); err_eye_tgt_py{counter_session,1}{counter_harvest,1}(1,counter_sac)];

                            data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1}(1,counter_sac) = err_rot(1);
                            data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1}(1,counter_sac) = err_rot(2);
                            data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1}(1,counter_sac) = sqrt(data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1}(1,counter_sac)^2 + data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1}(1,counter_sac)^2);
                            %  data.det_var_err_eye_tgt_{counter_session,1}{counter_harvest,1}(1,counter_sac) = det(nancov(data.err_eye_tgt_px_rot_{counter_session,1}{counter_harvest,1}(1,counter_sac), data.err_eye_tgt_py_rot_{counter_session,1}{counter_harvest,1}(1,counter_sac)));
                        end
                        if isempty(counter_sac)
                            data.err_eye_tgt_px_rot{counter_session,1}{counter_harvest,1} = [];
                            data.err_eye_tgt_py_rot{counter_session,1}{counter_harvest,1} = [];
                            data.dist_err_eye_tgt{counter_session,1}{counter_harvest,1} = [];
                        end
                    end
                end
            end
        end
    end
end

end
%% function pad_data
function [data] = pad_data(data, count)
fieldname = fieldnames(data);
for counter_session = 1 : count.num_session
    for counter_type = 1 : count.num_type
        if counter_type == 1
            for counter_fieldname = 1 : length(fieldname)
                % check if lick or sac related data
                if sum(contains(string(fieldname),'_lick')) > 0
                    % check if harvest related data or lick related data
                    if ~contains(string(fieldname(counter_fieldname)),[ "_harvest"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) < 500
                        for counter_harvest = 1 : count.num_harvest(counter_session,1)
                            pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest{counter_session,1}(counter_harvest,1);
                            if pad_val_within >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[0 pad_val_within], nan, 'post');
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_lick_harvest);
                            end
                        end
                        data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                    elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},2) == 500
                        for counter_harvest = 1 : count.num_harvest(counter_session,1)
                            pad_val_within = count.max_num_lick_harvest - count.num_lick_harvest{counter_session,1}(counter_harvest,1);
                            if pad_val_within >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[pad_val_within 0], nan, 'post');
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(1:1:count.max_num_lick_harvest,:);
                            end
                        end
                    elseif contains(string(fieldname(counter_fieldname)),["_harvest"]) && ~contains(string(fieldname(counter_fieldname)),["tag"])
                        pad_val_accross = count.max_num_harvest - count.num_harvest(counter_session,1);
                        if pad_val_accross >= 0
                            data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                        else
                            data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest,1))';
                        end
                    end

                elseif sum(contains(string(fieldname),'_sac')) > 0
                    % check if work related data or sac related data
                    if ~contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work", "num_trial_attempt_work"]) && size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) < 150
                        for counter_harvest = 1 : count.num_harvest(counter_session,1)
                            if(contains(string(fieldname(counter_fieldname)),["time_start_trial","time_end_trial", "ITI","ITR"]))
                                pad_val_within = count.max_num_trial_work - count.num_trial_work{counter_session,1}{counter_harvest,1};
                            elseif(contains(string(fieldname(counter_fieldname)),["_h"])) && ~(contains(string(fieldname(counter_fieldname)),["fix", "before", "after"]))
                                pad_val_within = count.max_num_sac_work - count.num_sac_harvest{counter_session,1}(counter_harvest,1);
                            elseif(contains(string(fieldname(counter_fieldname)),["_h", "fix", "before", "after"]))
                                pad_val_within = count.max_num_sac_work - count.num_fix_harvest{counter_session,1}(counter_harvest,1);
                            else
                                pad_val_within = count.max_num_sac_work - count.num_sac_work{counter_session,1}(counter_harvest,1);
                            end
                            if isempty(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1})
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = nan(1, pad_val_within);
                            else
                                if pad_val_within >= 0
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1},[0 pad_val_within], nan, 'post');
                                else
                                    data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_sac_work);
                                end
                            end
                        end
                        data.(string(fieldname(counter_fieldname))){counter_session,1} =cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1});
                    elseif size(data.(string(fieldname(counter_fieldname))){counter_session,1}{1,:},1) == 150
                        for counter_harvest = 1 : count.num_harvest(counter_session,1)
                            if(contains(string(fieldname(counter_fieldname)),["_h"]))
                                pad_val_within = count.max_num_sac_work - count.num_sac_harvest{counter_session,1}(counter_harvest,1);
                            else
                                pad_val_within = count.max_num_sac_work - count.num_sac_work{counter_session,1}(counter_harvest,1);
                            end
                            if pad_val_within >= 0
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = padarray(data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}',[pad_val_within 0], nan, 'post');
                            else
                                data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,:} = data.(string(fieldname(counter_fieldname))){counter_session,1}{counter_harvest,1}(:,1:1:count.max_num_sac_work)';
                            end
                        end
                    elseif contains(string(fieldname(counter_fieldname)),["duration_work", "num_sac_work", "rew_gained_work", "num_trial_work", "num_trial_attempt_work"])
                        pad_val_accross = count.max_num_harvest - count.num_harvest(counter_session,1);
                        if pad_val_accross >= 0
                            data.(string(fieldname(counter_fieldname))){counter_session,1} = padarray(cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}),[pad_val_accross 0], nan, 'post')';
                        else
                            data.(string(fieldname(counter_fieldname))){counter_session,1} = cell2mat(data.(string(fieldname(counter_fieldname))){counter_session,1}(1:1:count.max_num_harvest,1))';
                        end
                    end
                end
            end
        end
    end
end
end
%% function build bins
function [count] = build_bins(lick, sac, count, path)
%%% interval %%%
interval_bin_ = [0.2 0.4 0.6 0.8 1]; % percentages of total intervals

%%% rew %%%
rew_edges = [0 0.25 0.5 0.75 1 1.25]; % rew_bins = [0.25 0.5 0.75 1 1.25];

%%% trials %%%
num_trial_edges = [0 2 4 6 8 10]; % num trials completed in previous work interval

%%% rew rate %%%
%[~,rew_cosumed_rate_edges,~] = histcounts(cell2mat(lick.rew_consumed_rate_harvest),5)
rew_cosumed_rate_edges = [0 0.01 0.02 0.03 0.04 0.1 ];
rew_acquired_rate_edges = [0 0.01 0.02 0.03 0.04 0.1 ];

%%% weight %%%
weight_edges = [(min(count.weight)-1) median(count.weight) (max(count.weight)+1)];
%[~,weight_edges,~] = histcounts(count.weight,3);
% if contains(path.path_data_monkey_sorted,'59d') % animal weight before task
%     weight_edges = [330 345 355 370];
% elseif contains(path.path_data_monkey_sorted,'125d')
%     weight_edges = [320 340 350 370];
% end

%%% dist %%%
%[~,dist_edges,~] = histcounts(cell2mat(count.dist_tube_harvest_mean),11);
if contains(path.path_data_monkey_sorted,'59d')
    dist_edges = [7.2 9.78 10.64 11.93];
elseif contains(path.path_data_monkey_sorted,'125d')
    dist_edges = [7.6 9.9 10.82 12.66];
end
% dist_edges = [7 9.5 10.5 13]; % distances of tube edge centroid  from 0.
% dist_edges = [(min(cell2mat(count.dist_tube_harvest_mean))-1) median(cell2mat(count.dist_tube_harvest_mean)) (max(cell2mat(count.dist_tube_harvest_mean))-1)]; % distances of tube edge centroid  from 0.

% build
for counter_session = 1 : count.num_session
    % initialize
    interval_bin{counter_session,1} = nan(count.num_harvest(counter_session),1);
    num_interval_bin(counter_session,:) = floor(count.num_harvest(counter_session) * interval_bin_);

    % bin intervals
    for counter_bin = 1 : length(num_interval_bin(counter_session,:))
        if counter_bin == 1
            count.interval_bin_{counter_session,1}(1:num_interval_bin(counter_session,1),1) = counter_bin;
        else
            count.interval_bin_{counter_session,1}(num_interval_bin(counter_session,counter_bin-1) :num_interval_bin(counter_session,counter_bin),1) = counter_bin;
        end
    end

    % bin session day
    count.session_bin_{counter_session,1} = counter_session * ones(count.num_harvest(counter_session),1);

    % bin rew
    for counter_harvest = 1 : count.num_harvest(counter_session)
        count.rew_available_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(lick.rew_str_harvest{counter_session,1}(1:count.num_harvest(counter_session,1)), rew_edges);
        count.rew_consumed_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(lick.rew_consumed_harvest{counter_session,1}(1:count.num_harvest(counter_session,1)), rew_edges);
        count.rew_gained_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(sac.rew_gained_work{counter_session,1}(1:count.num_harvest(counter_session,1)), rew_edges);
        count.capture_rate_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(lick.rew_consumed_rate_harvest {counter_session,1}(1:count.num_harvest(counter_session,1)), rew_cosumed_rate_edges);
        count.acquire_rate_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(lick.rew_acquired_rate_harvest {counter_session,1}(1:count.num_harvest(counter_session,1)), rew_acquired_rate_edges);
    end
    count.rew_history_bin_{counter_session,1} = [nan;count.rew_consumed_bin_{counter_session,1}(1:end-1)];
    count.capture_rate_history_bin_{counter_session,1} = [nan;count.capture_rate_bin_{counter_session,1}(1:end-1)];
    count.acquire_rate_history_bin_{counter_session,1} = [nan;count.acquire_rate_bin_{counter_session,1}(1:end-1)];

    % bin dist tube
    %         for counter_harvest = 1 : count.num_harvest(counter_session)
    %             count.dist_tube_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(lick.dist_tube_harvest{counter_session,1}(1:count.num_harvest(counter_session,1)), dist_edges);
    %         end
    %     mode_dist_edges = mode(lick.dist_tube_harvest{counter_session, 1});
    %    count.dist_tube_bin_{counter_session,1} = discretize(count.dist_tube_harvest(counter_session), dist_edges) *ones(count.num_harvest(counter_session,1),1);
    count.dist_tube_bin_{counter_session,1} = discretize(count.dist_tube_harvest_mean{counter_session,1}, dist_edges) ;

    % bin num trial
    for counter_harvest = 1 : count.num_harvest(counter_session)
        count.num_trial_bin_{counter_session,1}(1:count.num_harvest(counter_session,1),1) = discretize(sac.num_trial_work{counter_session,1}(1:count.num_harvest(counter_session,1)), num_trial_edges);
    end

    % bin weight
    count.weight_bin_{counter_session,1} = discretize(count.weight(counter_session,1), weight_edges) *ones(count.num_harvest(counter_session,1),1);

    % bin harvest direction
    count.dir_harvest_bin_{counter_session,1}(count.dir_harvest{counter_session,1} == 'r',1) = 1;
    count.dir_harvest_bin_{counter_session,1}(count.dir_harvest{counter_session,1} == 'l',1) = 2;
    count.dir_harvest_bin_{counter_session,1}(count.dir_harvest{counter_session,1} == 'N',1) = 0;

    % pad bin data
    pad_val_accross = count.max_num_harvest - count.num_harvest(counter_session,1);
    count.interval_bin{counter_session,1} = padarray(count.interval_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.session_bin{counter_session,1} = padarray(count.session_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.rew_available_bin{counter_session,1} = padarray(count.rew_available_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.rew_consumed_bin{counter_session,1} = padarray(count.rew_consumed_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.rew_gained_bin{counter_session,1} = padarray(count.rew_gained_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.rew_history_bin{counter_session,1} = padarray(count.rew_history_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.dist_tube_bin{counter_session,1} = padarray(count.dist_tube_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.num_trial_bin{counter_session,1} = padarray(count.num_trial_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.weight_bin{counter_session,1} = padarray(count.weight_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.dir_harvest_bin{counter_session,1} = padarray(count.dir_harvest_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.capture_rate_bin{counter_session,1} = padarray(count.capture_rate_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.capture_rate_history_bin{counter_session,1} = padarray(count.capture_rate_history_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.acquire_rate_bin{counter_session,1} = padarray(count.acquire_rate_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
    count.acquire_rate_history_bin{counter_session,1} = padarray(count.acquire_rate_history_bin_{counter_session,1},[pad_val_accross 0], nan, 'post');
end
end

%% function add_RPE
function [count] = add_RPE(lick, count)
for counter_session = 1 : count.num_session
    % initialize
    for counter = 1 : 2 % for prev and current
        count.is_TT{counter_session,counter} = logical(zeros(count.num_harvest(counter_session),count.max_num_lick_harvest));
        count.is_FF{counter_session,counter} = logical(zeros(count.num_harvest(counter_session),count.max_num_lick_harvest));
        count.is_FT{counter_session,counter} = logical(zeros(count.num_harvest(counter_session),count.max_num_lick_harvest));
        count.is_TF{counter_session,counter} = logical(zeros(count.num_harvest(counter_session),count.max_num_lick_harvest));
    end
    for counter_harvest = 1 : count.num_harvest(counter_session)
        for counter_lick = 2 : count.max_num_lick_harvest
            % tag 3 as F
            %             count.is_TT{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 6) && ...
            %                 (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 6);
            %             count.is_FF{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 7) && ...
            %                 (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 7);
            %             count.is_FT{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 6) && ...
            %                 (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 7);
            %             count.is_TF{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 7) && ...
            %                 (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 6);

            % tag 3 as T
            count.is_TT{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 6) && ...
                (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 6);
            count.is_FF{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 7) && ...
                (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 7);
            count.is_FT{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 3|| lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 6) && ...
                (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 7);
            count.is_TF{counter_session,2}(counter_harvest,counter_lick) = (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 5 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick) == 7) && ...
                (lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 2 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 3 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 4 || lick.tag_lick{counter_session, 1}(counter_harvest,counter_lick - 1) == 6);
        end
        count.is_TT{counter_session,1}(counter_harvest, find(count.is_TT{counter_session,2}(counter_harvest,:) == 1) - 1) = 1;
        count.is_FF{counter_session,1}(counter_harvest, find(count.is_FF{counter_session,2}(counter_harvest,:) == 1) - 1) = 1;
        count.is_FT{counter_session,1}(counter_harvest, find(count.is_FT{counter_session,2}(counter_harvest,:) == 1) - 1) = 1;
        count.is_TF{counter_session,1}(counter_harvest, find(count.is_TF{counter_session,2}(counter_harvest,:) == 1) - 1) = 1;
    end
end
end
%% function compute_vigor
function [FIT] = compute_vigor(data_lick, data_sac, count, path, run)
if run == 1
    % build and concatenate amp and vel data
    for counter_session = 1 : count.num_session
        tag_lick{1,counter_session} = data_lick.tag{counter_session, 1};
        tongue_dm_max{1,counter_session} = data_lick.tongue_dm_max{counter_session, 1};
        tongue_vm_max{1,counter_session} = data_lick.tongue_vm_max{counter_session, 1};
        tongue_vm_min{1,counter_session} = abs(data_lick.tongue_vm_min{counter_session, 1});

        tag_sac{1,counter_session} = data_sac.tag{counter_session, 1};
        eye_dm_max{1,counter_session} = data_sac.eye_r_amp_m{counter_session, 1};
        eye_vm_max{1,counter_session} =  max(data_sac.eye_r_vm{counter_session, 1});
    end
    %     clearvars data_lick data_sac

    tongue_dm_max = cell2mat(tongue_dm_max);
    tongue_vm_max = cell2mat(tongue_vm_max);
    tongue_vm_min = cell2mat(tongue_vm_min);
    tag_lick = cell2mat(tag_lick);

    eye_dm_max = cell2mat(eye_dm_max);
    eye_vm_max = cell2mat(eye_vm_max);
    tag_sac = cell2mat(tag_sac);

    % separate based on tag
    tag_lick_use = 1; % grooming
    tongue_dm_max_groom = tongue_dm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_max_groom  = tongue_vm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_min_groom  = tongue_vm_min(ismember(tag_lick,tag_lick_use));
    % is_outlier = isoutlier(tongue_dm_max_groom, 'quartiles') | isoutlier(tongue_vm_max_groom, 'quartiles') | isoutlier(tongue_vm_min_groom, 'quartiles');
    % tongue_dm_max_groom(is_outlier) = []; tongue_vm_max_groom(is_outlier) = []; tongue_vm_min_groom(is_outlier) = [];

    tag_lick_use = 2:7; % reward
    tongue_dm_max_rew = tongue_dm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_max_rew  = tongue_vm_max(ismember(tag_lick,tag_lick_use));
    tongue_vm_min_rew  = tongue_vm_min(ismember(tag_lick,tag_lick_use));
    %is_outlier = isoutlier(tongue_dm_max_rew, 'quartiles') | isoutlier(tongue_vm_max_rew, 'quartiles') | isoutlier(tongue_vm_min_rew, 'quartiles');
    %tongue_dm_max_rew(is_outlier) = []; tongue_vm_max_rew(is_outlier) = []; tongue_vm_min_rew(is_outlier) = [];

    tag_sac_use = 10; % irrelevant
    eye_dm_max_irr = eye_dm_max(ismember(tag_sac,tag_sac_use));
    eye_vm_max_irr  = eye_vm_max(ismember(tag_sac,tag_sac_use));
    %is_outlier = isoutlier(eye_dm_max_irr, 'quartiles') | isoutlier(eye_vm_max_irr, 'quartiles');
    %eye_dm_max_irr(is_outlier) = []; eye_vm_max_irr(is_outlier) = [];

    tag_sac_use = [1 4 6]; % relevant
    eye_dm_max_rel= eye_dm_max(ismember(tag_sac,tag_sac_use));
    eye_vm_max_rel = eye_vm_max(ismember(tag_sac,tag_sac_use));
    % is_outlier = isoutlier(eye_dm_max_rel, 'quartiles') | isoutlier(eye_vm_max_rel, 'quartiles');
    % eye_dm_max_rel(is_outlier) = []; eye_vm_max_rel(is_outlier) = [];

    % amplitude based bins for vigor plot
    % build tongue
    tongue_bins = [5:2.5:22.5];
    tongue_edges = [2.5:2.5:25];

    bin_tongue_groom = discretize(tongue_dm_max_groom,tongue_edges);
    bin_tongue_rew = discretize(tongue_dm_max_rew,tongue_edges);

    for counter_bin = 1:length(tongue_bins)
        mean_tongue_vm_max_groom_bin(counter_bin) = nanmean(tongue_vm_max_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_vm_max_groom_bin(counter_bin) = nanstd(tongue_vm_max_groom(bin_tongue_groom == counter_bin));%/sqrt(sum(bin_tongue_groom == counter_bin)) ;
        mean_tongue_vm_min_groom_bin(counter_bin) = nanmean(tongue_vm_min_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_vm_min_groom_bin(counter_bin) = nanstd(tongue_vm_min_groom(bin_tongue_groom == counter_bin));%/sqrt(sum(bin_tongue_groom == counter_bin)) ;
        mean_tongue_diff_groom_bin(counter_bin) = nanmean(tongue_vm_min_groom(bin_tongue_groom == counter_bin) - tongue_vm_max_groom(bin_tongue_groom == counter_bin)) ;
        sem_tongue_diff_groom_bin(counter_bin) = nanstd(tongue_vm_min_groom(bin_tongue_groom == counter_bin) - tongue_vm_max_groom(bin_tongue_groom == counter_bin));%/sqrt(sum(bin_tongue_groom == counter_bin)) ;

        mean_tongue_vm_max_rew_bin(counter_bin) = nanmean(tongue_vm_max_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_vm_max_rew_bin(counter_bin) = nanstd(tongue_vm_max_rew(bin_tongue_rew == counter_bin));%/sqrt(sum(bin_tongue_rew == counter_bin)) ;
        mean_tongue_vm_min_rew_bin(counter_bin) = nanmean(tongue_vm_min_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_vm_min_rew_bin(counter_bin) = nanstd(tongue_vm_min_rew(bin_tongue_rew == counter_bin));%/sqrt(sum(bin_tongue_rew == counter_bin)) ;
        mean_tongue_diff_rew_bin(counter_bin) = nanmean(tongue_vm_min_rew(bin_tongue_rew == counter_bin) - tongue_vm_max_rew(bin_tongue_rew == counter_bin)) ;
        sem_tongue_diff_rew_bin(counter_bin) = nanstd(tongue_vm_min_rew(bin_tongue_rew == counter_bin) - tongue_vm_max_rew(bin_tongue_rew == counter_bin));%/sqrt(sum(bin_tongue_rew == counter_bin)) ;

        data_all_tongue_vm_max_{1, counter_bin} = tongue_vm_max_groom(bin_tongue_groom == counter_bin)';
        data_all_tongue_vm_max_{2, counter_bin} = tongue_vm_max_rew(bin_tongue_rew == counter_bin)';
        data_all_tongue_vm_min_{1, counter_bin} = tongue_vm_min_groom(bin_tongue_groom == counter_bin)';
        data_all_tongue_vm_min_{2, counter_bin} = tongue_vm_min_rew(bin_tongue_rew == counter_bin)';
        data_all_tongue_diff_{1, counter_bin} = (tongue_vm_min_groom(bin_tongue_groom == counter_bin) - tongue_vm_max_groom(bin_tongue_groom == counter_bin))';
        data_all_tongue_diff_{2, counter_bin} = (tongue_vm_min_rew(bin_tongue_rew == counter_bin) - tongue_vm_max_rew(bin_tongue_rew == counter_bin))';
        data_all_tongue_diff_groom_{1, counter_bin} = data_all_tongue_diff_{1, counter_bin};
        data_all_tongue_diff_rew_{1, counter_bin} = data_all_tongue_diff_{2, counter_bin};
    end

    % remove empty and corresponding bins
    position = 1;
    for counter_x = 1 : size(data_all_tongue_vm_max_,2)
        if ~isempty(data_all_tongue_vm_max_{1,counter_x}) && ~isempty(data_all_tongue_vm_max_{2,counter_x})
            data_all_tongue_vm_max{1,position} = data_all_tongue_vm_max_{1,counter_x}; data_all_tongue_vm_max{2,position} = data_all_tongue_vm_max_{2,counter_x};
            data_all_tongue_vm_min{1,position} = data_all_tongue_vm_min_{1,counter_x}; data_all_tongue_vm_min{2,position} = data_all_tongue_vm_min_{2,counter_x};
            data_all_tongue_diff{1,position} = data_all_tongue_diff_{1,counter_x}; data_all_tongue_diff{2,position} = data_all_tongue_diff_{2,counter_x};
            data_all_tongue_diff_groom{1,position} = data_all_tongue_diff_groom_{1,counter_x};
            data_all_tongue_diff_rew{1,position} = data_all_tongue_diff_rew_{1,counter_x};

            position = position + 1;
        end
    end

    [test_tongue_vm_max, summary_tongue_vm_max,table_data] = PGH_stat_ANOVA(data_all_tongue_vm_max);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_tongue_vm_max' '.csv']);

    [test_tongue_vm_min, summary_tongue_vm_min,table_data] = PGH_stat_ANOVA(data_all_tongue_vm_min);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_tongue_vm_min' '.csv']);

    [test_tongue_diff, summary_tongue_diff,table_data] = PGH_stat_ANOVA(data_all_tongue_diff);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_tongue_diff' '.csv']);

    [test_tongue_diff_groom, summary_tongue_diff_groom,table_data] = PGH_stat_ANOVA(data_all_tongue_diff_groom);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_tongue_diff_groom' '.csv']);

    [test_tongue_diff_rew, summary_tongue_diff_rew,table_data] = PGH_stat_ANOVA(data_all_tongue_diff_rew);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_tongue_diff_rew' '.csv']);


    % build eye
    eye_bins = [1:10];
    eye_edges = [0:11];

    bin_eye_irr = discretize(eye_dm_max_irr,eye_edges);
    bin_eye_rel = discretize(eye_dm_max_rel,eye_edges);

    for counter_bin = 1:length(eye_bins)
        mean_eye_vm_max_irr_bin(counter_bin) = nanmean(eye_vm_max_irr(bin_eye_irr == counter_bin)) ;
        sem_eye_vm_max_irr_bin(counter_bin) = nanstd(eye_vm_max_irr(bin_eye_irr == counter_bin));%/sqrt(sum(bin_eye_irr == counter_bin)) ;

        mean_eye_vm_max_rel_bin(counter_bin) = nanmean(eye_vm_max_rel(bin_eye_rel == counter_bin)) ;
        sem_eye_vm_max_rel_bin(counter_bin) = nanstd(eye_vm_max_rel(bin_eye_rel == counter_bin));%/sqrt(sum(eye_vm_max_rel == counter_bin)) ;

        data_all_eye_vm_max_{1, counter_bin} = eye_vm_max_irr(bin_eye_irr == counter_bin)';
        data_all_eye_vm_max_{2, counter_bin} = eye_vm_max_rel(bin_eye_rel == counter_bin)';
    end

    % remove empty and corresponding bins
    position = 1;
    for counter_x = 1 : size(data_all_eye_vm_max_,2)
        if ~isempty(data_all_eye_vm_max_{1,counter_x}) && ~isempty(data_all_eye_vm_max_{2,counter_x})
            data_all_eye_vm_max{1,position} = data_all_eye_vm_max_{1,counter_x}; data_all_eye_vm_max{2,position} = data_all_eye_vm_max_{2,counter_x};
            position = position + 1;
        end
    end

    [test_eye_vm_max, summary_eye_vm_max, table_data] = PGH_stat_ANOVA(data_all_eye_vm_max);
    writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted 'vigor_eye_vm_max'  '.csv' ]);


    % build fits
    B0 = [100; 0.1];
    f = fittype('a*(1-(1./(1+b*x)))');

    FIT.fit_lick_groom_pro = fit(tongue_dm_max_groom',tongue_vm_max_groom',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_lick_rew_pro = fit(tongue_dm_max_rew',tongue_vm_max_rew',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_lick_pro = fit([tongue_dm_max_groom tongue_dm_max_rew]',[tongue_vm_max_groom tongue_vm_max_rew]',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);

    FIT.fit_lick_groom_ret = fit(tongue_dm_max_groom',tongue_vm_min_groom',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_lick_rew_ret = fit(tongue_dm_max_rew',tongue_vm_min_rew',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_lick_ret = fit([tongue_dm_max_groom tongue_dm_max_rew]',[tongue_vm_min_groom tongue_vm_min_rew]',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);

    FIT.fit_sac_irr = fit(eye_dm_max_irr',eye_vm_max_irr',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_sac_rel = fit(eye_dm_max_rel',eye_vm_max_rel',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);
    FIT.fit_eye = fit([eye_dm_max_irr eye_dm_max_rel]',[eye_vm_max_irr eye_vm_max_rel]',f,'StartPoint',B0,'Robust','on','MaxIter',1e3,'Lower',[0,0]);

    save([path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_FIT.mat' ], 'FIT');

    % plot vigor
    plot_vigor = 1;
    if plot_vigor == 1
        % plot
        fig_vigor = figure;
        s1 = subplot(3,6,1); %%% FIT LICK PRO
        hold on
        plot(tongue_dm_max_groom,tongue_vm_max_groom,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_groom_pro, 'r')
        plot(tongue_dm_max_rew,tongue_vm_max_rew,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_rew_pro, 'b')
        plot(FIT.fit_lick_pro, 'k')
        xlabel('lick amp. (mm)')
        ylabel('lick pro. speed (mm/s)')
        legend('off')
        title('groom (r) vs rew (b)')
        ylim([0 600])
        yticks([0 : 100: 600])
        xlim([0 25])
        xticks([0 : 5: 25])

        s2 = subplot(3,6,2); %%% STAT LICK PRO
        hold on
        tongue_vigor_pro_groom = tongue_vm_max_groom'./FIT.fit_lick_pro(tongue_dm_max_groom);
        tongue_vigor_ret_groom = tongue_vm_min_groom'./FIT.fit_lick_ret(tongue_dm_max_groom);
        tongue_vigor_pro_rew = tongue_vm_max_rew'./FIT.fit_lick_pro(tongue_dm_max_rew);
        tongue_vigor_ret_rew = tongue_vm_min_rew'./FIT.fit_lick_ret(tongue_dm_max_rew);

        tongue_vigor_pro_groom(isoutlier(tongue_vigor_pro_groom, 'quartiles') | isoutlier(tongue_vigor_ret_groom, 'quartiles')) = nan;
        tongue_vigor_pro_rew(isoutlier(tongue_vigor_pro_rew, 'quartiles') | isoutlier(tongue_vigor_ret_rew, 'quartiles')) = nan;
        % t test
        [h,p,ci,stats] = ttest2(tongue_vigor_pro_groom, tongue_vigor_pro_rew, "Tail", "left");
        % 95CI SE*1.96
        mean_tongue_vigor_pro_groom = nanmean(tongue_vigor_pro_groom);
        CI_tongue_vigor_pro_groom = nanstd(tongue_vigor_pro_groom)/sqrt(length(tongue_vigor_pro_groom))*1.96;
        mean_tongue_vigor_pro_rew = nanmean(tongue_vigor_pro_rew);
        CI_tongue_vigor_pro_rew = nanstd(tongue_vigor_pro_rew)/sqrt(length(tongue_vigor_pro_rew))*1.96;
        test = ['2s/1t t-test, t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', '  num2str(p) newline...
            'mean+/-95%CI, groom: ' num2str(mean_tongue_vigor_pro_groom) '+/-' num2str(CI_tongue_vigor_pro_groom) newline ...
            ', rew: ' num2str(mean_tongue_vigor_pro_rew) '+/-' num2str(CI_tongue_vigor_pro_rew)];
        x_axis = {'groom', 'rew'};
        grp = [zeros(1,length(tongue_vigor_pro_groom)), ones(1,length(tongue_vigor_pro_rew))]';
        %         boxplot([tongue_vigor_pro_groom; tongue_vigor_pro_rew],grp,'Labels',x_axis,'Colors','k','whisker', inf);
        %         errorbar([1 2],[mean_tongue_vigor_pro_groom mean_tongue_vigor_pro_rew], [CI_tongue_vigor_pro_groom CI_tongue_vigor_pro_rew], '.r', 'LineWidth', 2, 'MarkerSize', 10);
        x_labels = {}; x_labels(grp == 0) = x_axis(1); x_labels(grp == 1) = x_axis(2);
        violinplot([tongue_vigor_pro_groom; tongue_vigor_pro_rew],x_labels, 'ShowMedian', false, 'ShowMean',true,'ShowData', false);
        ylabel('lick pro. vigor')
        %title(test, 'Interpreter', 'none')

        s3 = subplot(3,6,3); %%% FIT LICK RET
        hold on
        plot(tongue_dm_max_groom,tongue_vm_min_groom,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_groom_ret, 'r')
        plot(tongue_dm_max_rew,tongue_vm_min_rew,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_lick_rew_ret, 'b')
        plot(FIT.fit_lick_ret, 'k')

        xlabel('lick amp. (mm)')
        ylabel('lick ret. speed (mm/s)')
        legend('off')
        title('groom (r) vs rew (b)')
        ylim([0 600])
        yticks([0 : 100: 600])
        xlim([0 25])
        xticks([0 : 5: 25])

        s4 = subplot(3,6,4); %%% STAT LICK RET
        hold on
        tongue_vigor_pro_groom = tongue_vm_max_groom'./FIT.fit_lick_pro(tongue_dm_max_groom);
        tongue_vigor_ret_groom = tongue_vm_min_groom'./FIT.fit_lick_ret(tongue_dm_max_groom);
        tongue_vigor_pro_rew = tongue_vm_max_rew'./FIT.fit_lick_pro(tongue_dm_max_rew);
        tongue_vigor_ret_rew = tongue_vm_min_rew'./FIT.fit_lick_ret(tongue_dm_max_rew);

        tongue_vigor_ret_groom(isoutlier(tongue_vigor_pro_groom, 'quartiles') | isoutlier(tongue_vigor_ret_groom, 'quartiles')) = nan;
        tongue_vigor_ret_rew(isoutlier(tongue_vigor_pro_rew, 'quartiles') | isoutlier(tongue_vigor_ret_rew, 'quartiles')) = nan;
        % t test
        [h,p,ci,stats] = ttest2(tongue_vigor_ret_groom, tongue_vigor_ret_rew, "Tail", "left");
        % 95CI SE*1.96
        mean_tongue_vigor_ret_groom = nanmean(tongue_vigor_ret_groom);
        CI_tongue_vigor_ret_groom = nanstd(tongue_vigor_ret_groom)/sqrt(length(tongue_vigor_ret_groom))*1.96;
        mean_tongue_vigor_ret_rew = nanmean(tongue_vigor_ret_rew);
        CI_tongue_vigor_ret_rew = nanstd(tongue_vigor_ret_rew)/sqrt(length(tongue_vigor_ret_rew))*1.96;
        test = ['2s/1t t-test, t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', '  num2str(p) newline...
            'mean+/-95%CI, groom: ' num2str(mean_tongue_vigor_ret_groom) '+/-' num2str(CI_tongue_vigor_ret_groom) newline ...
            ', rew: ' num2str(mean_tongue_vigor_ret_rew) '+/-' num2str(CI_tongue_vigor_ret_rew)];
        x_axis = {'groom', 'rew'};
        grp = [zeros(1,length(tongue_vigor_ret_groom)), ones(1,length(tongue_vigor_ret_rew))]';
        %         boxplot([tongue_vigor_ret_groom; tongue_vigor_ret_rew],grp,'Labels',x_axis,'Colors','k','whisker', inf);
        %         errorbar([1 2],[mean_tongue_vigor_ret_groom mean_tongue_vigor_ret_rew], [CI_tongue_vigor_ret_groom CI_tongue_vigor_ret_rew], '.r', 'LineWidth', 2, 'MarkerSize', 10);
        x_labels = {}; x_labels(grp == 0) = x_axis(1); x_labels(grp == 1) = x_axis(2);
        violinplot([tongue_vigor_ret_groom; tongue_vigor_ret_rew],x_labels, 'ShowMedian', false, 'ShowMean',true,'ShowData', false);
        ylabel('lick ret. vigor')
        %title(test, 'Interpreter', 'none')

        s5 = subplot(3,6,5); %%% FIT SAC
        hold on
        plot(eye_dm_max_irr,eye_vm_max_irr,'or', 'MarkerSize', 0.1)
        plot(FIT.fit_sac_irr, 'r')
        plot(eye_dm_max_rel,eye_vm_max_rel,'ob', 'MarkerSize', 0.1)
        plot(FIT.fit_sac_rel, 'b')
        plot(FIT.fit_eye, 'k')
        xlabel('saccade amp. (deg)')
        ylabel('saccade vel. (deg/s)')
        legend('off')
        title('task irr (r) vs task rel (b)')

        s6 = subplot(3,6,6); %%% STAT SAC
        hold on
        eye_vigor_irr = eye_vm_max_irr'./FIT.fit_eye(eye_dm_max_irr);eye_vigor_irr(isoutlier(eye_vigor_irr, 'quartiles')) = nan;
        eye_vigor_rel = eye_vm_max_rel'./FIT.fit_eye(eye_dm_max_rel);eye_vigor_rel(isoutlier(eye_vigor_rel, 'quartiles')) = nan;
        % t test
        [h,p,ci,stats] = ttest2(eye_vigor_irr, eye_vigor_rel, "Tail", "left");
        % 95CI SE*1.96
        mean_eye_vigor_irr = nanmean(eye_vigor_irr);
        CI_eye_vigor_irr = nanstd(eye_vigor_irr)/sqrt(length(eye_vigor_irr))*1.96;
        mean_eye_vigor_rel = nanmean(eye_vigor_rel);
        CI_eye_vigor_rel = nanstd(eye_vigor_rel)/sqrt(length(eye_vigor_rel))*1.96;
        test = ['2s/1t t-test, t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', '  num2str(p) newline...
            'mean+/-95%CI, irr: ' num2str(mean_eye_vigor_irr) '+/-' num2str(CI_eye_vigor_irr) newline ...
            ', rel: ' num2str(mean_eye_vigor_rel) '+/-' num2str(CI_eye_vigor_rel)];
        x_axis = {'irr', 'rel'};
        grp = [zeros(1,length(eye_vigor_irr)), ones(1,length(eye_vigor_rel))]';
        % boxplot([eye_vigor_irr; eye_vigor_rel],grp,'Labels',x_axis,'Colors','k','whisker', inf);
        % errorbar([1 2],[mean_eye_vigor_irr mean_eye_vigor_rel], [CI_eye_vigor_irr CI_eye_vigor_rel], '.r', 'LineWidth', 2, 'MarkerSize', 10);
        x_labels = {}; x_labels(grp == 0) = x_axis(1); x_labels(grp == 1) = x_axis(2);
        violinplot([eye_vigor_irr; eye_vigor_rel],x_labels, 'ShowMedian', false, 'ShowMean',true, 'ShowData', false);
        ylabel('saccade vigor')
        %title(test, 'Interpreter', 'none')

        s7_8 = subplot(3,6,[7 8]); %%% PLOT LICK PRO VIG
        hold on
        plot(tongue_bins,mean_tongue_vm_max_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_vm_max_groom_bin, sem_tongue_vm_max_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_vm_max_groom_bin + sem_tongue_vm_max_groom_bin, 'r', tongue_bins, mean_tongue_vm_max_groom_bin - sem_tongue_vm_max_groom_bin, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(tongue_bins,mean_tongue_vm_max_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_vm_max_rew_bin, sem_tongue_vm_max_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_vm_max_rew_bin + sem_tongue_vm_max_rew_bin, 'b', tongue_bins, mean_tongue_vm_max_rew_bin - sem_tongue_vm_max_rew_bin, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        xlabel('lick amp. (mm)')
        ylabel('lick pro. speed (mm/s)')
        title(test_tongue_vm_max, 'Interpreter', 'none')

        s9_10 = subplot(3,6,[9 10]); %%% PLOT LICK RET VIG
        hold on
        plot(tongue_bins,mean_tongue_vm_min_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_vm_min_groom_bin, sem_tongue_vm_min_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_vm_min_groom_bin + sem_tongue_vm_min_groom_bin, 'r', tongue_bins, mean_tongue_vm_min_groom_bin - sem_tongue_vm_min_groom_bin, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(tongue_bins,mean_tongue_vm_min_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_vm_min_rew_bin, sem_tongue_vm_min_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_vm_min_rew_bin + sem_tongue_vm_min_rew_bin, 'b', tongue_bins, mean_tongue_vm_min_rew_bin - sem_tongue_vm_min_rew_bin, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        xlabel('lick amp. (mm)')
        ylabel('lick ret. speed (mm/s)')
        ylim([0 600])
        yticks([0 : 100: 600])
        xlim([0 25])
        xticks([0 : 5: 25])
        title(test_tongue_vm_min, 'Interpreter', 'none')

        s11_12 = subplot(3,6,[11 12]); %% PLOT SAC VIG
        hold on
        plot(eye_bins,mean_eye_vm_max_irr_bin,'r')
        %shade(eye_bins, mean_eye_vm_max_irr_bin + sem_eye_vm_max_irr_bin, 'r', eye_bins, mean_eye_vm_max_irr_bin - sem_eye_vm_max_irr_bin, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        errorbar(eye_bins,mean_eye_vm_max_irr_bin, sem_eye_vm_max_irr_bin, '.r', 'LineWidth', 2, 'MarkerSize', 10);
        plot(eye_bins,mean_eye_vm_max_rel_bin,'b')
        errorbar(eye_bins,mean_eye_vm_max_rel_bin, sem_eye_vm_max_rel_bin, '.b', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(eye_bins, mean_eye_vm_max_rel_bin + sem_eye_vm_max_rel_bin, 'b', eye_bins, mean_eye_vm_max_rel_bin - sem_eye_vm_max_rel_bin, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        xlabel('lick amp. (mm)')
        xlabel('saccade amp. (deg)')
        ylabel('saccade vel. (deg/s)')
        title(test_eye_vm_max, 'Interpreter', 'none')

        s13_14 = subplot(3,6,[13 14]); %%% PLOT LICK DIFF
        hold on
        plot(tongue_bins,mean_tongue_diff_groom_bin,'r')
        errorbar(tongue_bins,mean_tongue_diff_groom_bin, sem_tongue_diff_groom_bin, '.r', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_diff_groom_bin + sem_tongue_diff_groom_bin, 'r', tongue_bins, mean_tongue_diff_groom_bin - sem_tongue_diff_groom_bin, 'r', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','r', 'FillColor', 'r');
        plot(tongue_bins,mean_tongue_diff_rew_bin,'b')
        errorbar(tongue_bins,mean_tongue_diff_rew_bin, sem_tongue_diff_rew_bin, '.b', 'LineWidth', 2, 'MarkerSize', 10);
        %shade(tongue_bins, mean_tongue_diff_rew_bin + sem_tongue_diff_rew_bin, 'b', tongue_bins, mean_tongue_diff_rew_bin - sem_tongue_diff_rew_bin, 'b', 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color','b', 'FillColor', 'b');
        xlabel('lick amp. (mm)')
        ylabel('diff. ret.-pro. (mm/s)')
        ylim([-100 300])
        yticks([-100 : 50: 300])
        xlim([0 25])
        xticks([0 : 5: 25])
        title(test_tongue_diff, 'Interpreter', 'none')

        s15 = subplot(3,6,15); %%% STATS LICK DIFF
        hold on
        tongue_diff_groom = (tongue_vm_min_groom - tongue_vm_max_groom)';
        tongue_diff_rew = (tongue_vm_min_rew - tongue_vm_max_rew)';
        tongue_diff_groom(isoutlier(tongue_diff_groom, 'quartiles') | isoutlier(tongue_vm_min_groom', 'quartiles') | isoutlier(tongue_vm_max_groom', 'quartiles')) = nan;
        tongue_diff_rew(isoutlier(tongue_diff_rew, 'quartiles') | isoutlier(tongue_vm_min_rew', 'quartiles') | isoutlier(tongue_vm_max_rew', 'quartiles')) = nan;

        % t test
        [h_groom,p_groom,ci_groom,stats_groom] = ttest(tongue_diff_groom,0, "Tail", "right");
        [h_rew,p_rew,ci_rew,stats_rew] = ttest(tongue_diff_rew,0, "Tail", "right");
        % 95CI SE*1.96
        mean_tongue_diff_groom = nanmean(tongue_diff_groom);
        CI_tongue_diff_groom = nanstd(tongue_diff_groom)/sqrt(length(tongue_diff_groom))*1.96;
        mean_tongue_diff_rew = nanmean(tongue_diff_rew);
        CI_tongue_diff_rew = nanstd(tongue_diff_rew)/sqrt(length(tongue_diff_rew))*1.96;
        test = ['1s/1t t-test, groom: t(' num2str(stats_groom.df) ') = ' num2str(stats_groom.tstat) ', '  num2str(p_groom) ...
            ', mean+/-95%CI: ' num2str(mean_tongue_diff_groom) '+/-' num2str(CI_tongue_diff_groom) newline ...
            'rew: t(' num2str(stats_rew.df) ') = ' num2str(stats_rew.tstat) ', '  num2str(p_rew) ...
            ', mean+/-95%CI: ' num2str(mean_tongue_diff_rew) '+/-' num2str(CI_tongue_diff_rew) ];
        x_axis = {'groom', 'rew'};
        grp = [zeros(1,length(tongue_diff_groom)), ones(1,length(tongue_diff_rew))]';
        %         boxplot([tongue_diff_groom; tongue_diff_rew],grp,'Labels',x_axis,'Colors','k','whisker', inf);
        %         errorbar([1 2],[mean_tongue_diff_groom mean_tongue_diff_rew], [CI_tongue_diff_groom CI_tongue_diff_rew], '.r', 'LineWidth', 2, 'MarkerSize', 10);
        x_labels = {}; x_labels(grp == 0) = x_axis(1); x_labels(grp == 1) = x_axis(2);
        violinplot([tongue_diff_groom; tongue_diff_rew],x_labels, 'ShowMedian', false, 'ShowMean',true,'ShowData', false);
        ylabel('ret. - pro. speed (mm/s)')
        title(test, 'Interpreter', 'none')

        s16 = subplot(3,6,16); %%% STATS LICK DIFF
        title([test_tongue_diff_groom newline test_tongue_diff_rew], 'Interpreter', 'none')

        linkaxes([s1 s3 s7_8 s9_10], ['x','y']);
        linkaxes([s5 s11_12], ['x' ,'y']);

        sgtitle([path.path_data_monkey_sorted(6:end) ' | sessions: ' num2str(count.num_session) ' | licks (groom, rew): ' num2str(length(bin_tongue_groom)), ', ' num2str(length(bin_tongue_rew)) ' | sacs (irr, rel): ' num2str(length(bin_eye_irr)), ', ' num2str(length(bin_eye_rel))], 'Interpreter', 'none')
        ESN_Beautify_Plot
        fig_vigor.WindowState = 'maximized';
        saveas(gcf, [path.out_path 'VIGOR' filesep path.path_data_monkey_sorted '_vigor' ], 'pdf');
        %close(gcf)
    end

else
    FIT = [];
end
end
%% function filter_pupil
function [pupil_area] = filter_pupil(pupil_area_raw,counter_session,path)
pupil_area = pupil_area_raw;

% cap negative values at 0
pupil_area(pupil_area<0) = 0;

% find prom threshold using kmeans on all peaks
[peaks,~,~,~] = findpeaks(-pupil_area, 'Annotate', 'extent');
[~, cent] = kmeans(peaks', 2,"Start",[median(peaks) ; 0]);
threshold = abs(round((cent(1) + cent(2))/2));

% detect troughs
[~,ind_peak_neg,half_width,~] = findpeaks(-pupil_area, 'MinPeakHeight', -threshold,'Annotate', 'extent');
width = half_width; width(isoutlier(width)) = nanmean(width(~isoutlier(width))); width = round(width);
for counter_peak = 1 : length(ind_peak_neg)
    deletion_window = -width(counter_peak):1:width(counter_peak);
    inds_nan = ind_peak_neg(counter_peak) + deletion_window;
    inds_nan(inds_nan<=0 | inds_nan > length(pupil_area)) = [];
    pupil_area(inds_nan) = nan;
end
pupil_area(pupil_area < threshold) = nan;

plot_check = 1;
if plot_check == 1
    figure
    hold on
    plot(pupil_area_raw,'r')
    plot(pupil_area,'k')
    ylabel('pupil area')
    xlabel('time')
    ylim([0 inf])
    title(['session:' num2str(counter_session)])
    legend('raw', 'filtered')
    ESN_Beautify_Plot(gcf,[20,8])
    saveas(gcf, [path.out_path 'PUPIL' filesep path.path_data_monkey_sorted '_pupil_' num2str(counter_session) ], 'pdf');
    close all
end
end
%% function build_pupil
function [pupil_area_work, pupil_velocity_work, pupil_area_harvest, pupil_velocity_harvest] = build_pupil(pupil_area_work_, pupil_area_harvest_,time_1K_work_,time_1K_harvest_,time_lick_onset_,time_sac_onset_,count)
for counter_session = 1 : count.num_session
    for counter_harvest = 1 : count.num_harvest(counter_session)
        %pupil_area_harvest_{counter_session, 1}{counter_harvest, 1}(isoutlier(pupil_area_harvest_{counter_session, 1}{counter_harvest, 1})) = nan; % filter blinks and other artifacts
        %pupil_area_work_{counter_session, 1}{counter_harvest, 1}(isoutlier(pupil_area_work_{counter_session, 1}{counter_harvest, 1})) = nan; % filter blinks and other artifacts
        time_1K_harvest = time_1K_harvest_{counter_session, 1}{counter_harvest, 1};
        time_1K_work = time_1K_work_{counter_session, 1}{counter_harvest, 1};
        for counter_lick_harvest = 1 : count.num_lick_harvest{counter_session, 1}(counter_harvest, 1)
            time_lick_onset = time_lick_onset_{counter_session, 1}{counter_harvest, 1}(counter_lick_harvest);
            [~,ind] = min(abs(time_lick_onset-time_1K_harvest));
            if isempty(ind)
                pupil_area_harvest{counter_session, 1}{counter_harvest, 1}(counter_lick_harvest) = nan;
                pupil_velocity_harvest{counter_session, 1}{counter_harvest, 1}(counter_lick_harvest) = nan;
            else
                window = [-250:1:250] + ind; window(window<=0 | window>length(time_1K_harvest)) = []; % 500ms around lick
                pupil_area_harvest{counter_session, 1}{counter_harvest, 1}(counter_lick_harvest) = nanmean(pupil_area_harvest_{counter_session, 1}{counter_harvest, 1}(window));
                pupil_velocity_harvest{counter_session, 1}{counter_harvest, 1}(counter_lick_harvest) = nanmean(diff(pupil_area_harvest_{counter_session, 1}{counter_harvest, 1}(window)));
            end
        end

        for counter_sac_work = 1 : count.num_sac_work{counter_session, 1}(counter_harvest, 1)
            time_sac_onset = time_sac_onset_{counter_session, 1}{counter_harvest, 1}(counter_sac_work);
            [~,ind] = min(abs(time_sac_onset-time_1K_work));
            if isempty(ind)
                pupil_area_work{counter_session, 1}{counter_harvest, 1}(counter_sac_work) = nan;
                pupil_velocity_work{counter_session, 1}{counter_harvest, 1}(counter_sac_work) = nan;
            else
                window = [-250:1:250] + ind; window(window<=0 | window>length(time_1K_work)) = []; % 500ms around sac
                pupil_area_work{counter_session, 1}{counter_harvest, 1}(counter_sac_work) = nanmean(pupil_area_work_{counter_session, 1}{counter_harvest, 1}(window));
                pupil_velocity_work{counter_session, 1}{counter_harvest, 1}(counter_sac_work) = nanmean(diff(pupil_area_work_{counter_session, 1}{counter_harvest, 1}(window)));
            end
        end
        if count.num_sac_work{counter_session, 1}(counter_harvest, 1) == 0
            pupil_area_work{counter_session, 1}{counter_harvest, 1} = [];
            pupil_velocity_work{counter_session, 1}{counter_harvest, 1} = [];
        end
    end
end
end

%% function endpoint_error_var
function [bmean, bstd, data] = endpoint_error_var(data_x, data_y)
for counter_event = 1 : size(data_x,2)
    %data(:,counter_event) = det(nancov(data_x(:,counter_event), data_y(:,counter_event))) * ones(size(data_x,1),1);
    data(:,counter_event) = bootstrp(100,@(x)det(nancov(x)),[data_x(:,counter_event) data_y(:,counter_event)]);
    bmean(1,counter_event) = nanmean(data(:,counter_event));
    bstd(1,counter_event) = nanstd(data(:,counter_event));
end
end
%% function disclude_data
function [bin] = disclude_data(bin, disclude_l, disclude_r, dir_harvest_bin)
% left (dir == 2)
for counter_disclude_l = 1 : length(disclude_l)
    is_nan = dir_harvest_bin{disclude_l(counter_disclude_l), 1} == 2;
    bin{disclude_l(counter_disclude_l),1}(is_nan) = nan;
end
% right (dir == 1)
for counter_disclude_r = 1 : length(disclude_r)
    is_nan = dir_harvest_bin{disclude_r(counter_disclude_r), 1} == 1;
    bin{disclude_r(counter_disclude_r),1}(is_nan) = nan;
end
end
%% function PGH_stat_ANOVA
function [test,summary,table_data] = PGH_stat_ANOVA(data)
if size(data,1) == 1
    % bootstrap data
%     for counter_bin = 1 : length(data)
%         data_bs{1,counter_bin} = bootstrp(100,@(x)nanmean(x),data{1,counter_bin});
%     end
%     data = data_bs;

    % build data
    for counter_bin = 1 : length(data)
        length_data(counter_bin) = length(data{counter_bin});
    end
    max_length_data = max(length_data);
    t_data = nan(max_length_data, length(data));

    for counter_bin = 1 : length(data)
        t_data(1:length(data{1,counter_bin}),counter_bin) = data{1,counter_bin};
    end

    % for OSF
    table_data = array2table(t_data);

    [p, table_, stats] = anova1(t_data, [], 'off');
    [~,m,~,~] = multcompare(stats,'Display','off');
    [p_kw, table_kw, stats_kw] = kruskalwallis(t_data, [], 'off');
    test = ['anova1: F (' num2str(table_{2, 3}), ', ' num2str(table_{3, 3}), ') = ' num2str(round(table_{2, 5}, 2)) ', ' num2str(round(p,3,'significant')) newline ...
        'kw: X^2 (' num2str(table_kw{2, 3}), ', ' num2str(table_kw{3, 3}), ') = ' num2str(round(table_kw{2, 5}, 2)) ', ' num2str(round(p_kw,3,'significant'))];
    summary = round(m, 3, 'significant');

    if size(t_data,2) == 2
        [h, p, ~, stats] = ttest2(t_data(:,1), t_data(:,2));
        [p_wrs, h_wrs, stats_wrs] = ranksum(t_data(:,1), t_data(:,2));

        test = ['t(' num2str(stats.df) ') = ' num2str(round(stats.tstat, 3, 'significant')) ', p = '  num2str(round(p, 3,  'significant')) ', ' num2str(round(stats.sd, 3, 'significant')) newline ...
            'W = ' num2str(round(stats_wrs.zval, 3, 'significant')) ', p = '  num2str(round(p_wrs, 3,  'significant')) ', ' num2str(round(stats_wrs.ranksum, 3, 'significant'))];
    end
else
    % bootstrap data
%     for counter_bin = 1 : size(data, 1)
%         for counter_x = 1 : size(data,2)
%             data_bs{counter_bin,counter_x} = bootstrp(100,@(x)nanmean(x),data{counter_bin,counter_x});
%         end
%     end
%     data = data_bs;

    % build data
    for counter_bin = 1 : size(data, 1)
        for counter_x = 1 : size(data,2)
            length_data(counter_bin, counter_x) = length(data{counter_bin, counter_x});
        end
        max_length_data(counter_bin,1) =  max(length_data(counter_bin, :));

        for counter_x = 1 : size(data,2)
            t_data{counter_bin, counter_x} = nan(max_length_data(counter_bin,1),1);

            bin_id{counter_bin,counter_x} = counter_bin *ones(max_length_data(counter_bin,1),1);
            bin_x{counter_bin,counter_x} = counter_x *ones(max_length_data(counter_bin,1),1);
            t_data{counter_bin, counter_x}(1:length(data{counter_bin, counter_x})) = data{counter_bin, counter_x};
        end
    end
    bin_id = cell2mat(bin_id); bin_id = reshape(bin_id', numel(bin_id),1);
    bin_x = cell2mat(bin_x); bin_x = reshape(bin_x', numel(bin_x),1);
    t_data = cell2mat(t_data);t_data = reshape(t_data', numel(t_data),1);


    % for OSF
    table_data = array2table([bin_id, bin_x, t_data], 'VariableNames', {'bin_id', 'bin_x', 'data'});


    [p, table_, stats] = anovan(t_data,{bin_id bin_x},'model','full','varnames',{'bin_id','bin_x'}, 'display', 'off');
    [~,m,~,~] = multcompare(stats,'Display','off');
    test = ['bin_id: F (' num2str(table_{2, 3}), ', ' num2str(table_{5, 3}), ') = ' num2str(round(table_{2, 6}, 2)) ', ' num2str(round(table_{2, 7},3,'significant')) newline ...
        'bin_x: F (' num2str(table_{3, 3}), ', ' num2str(table_{5, 3}), ') = ' num2str(round(table_{3, 6}, 2)) ', ' num2str(round(table_{3, 7},3,'significant')) newline ...
        'bin_id*bin_x: F (' num2str(table_{4, 3}), ') = ' num2str(round(table_{4, 6}, 2)) ', ' num2str(round(table_{4, 7},3,'significant'))];
    summary = round(m, 3, 'significant');
end

end
%% function PGH_stat_RANOVA
function [test,summary, table_data] = PGH_stat_RANOVA(data)
% bootstrap data
% for counter_bin_dist = 1 : size(data, 1)
%     for counter_bin = 1 : size(data,2)
%         data_bs{counter_bin_dist,counter_bin} = bootstrp(100,@(x)nanmean(x),data{counter_bin_dist,counter_bin});
%     end
% end
% data = data_bs;

% build data
for counter_bin_dist = 1 : size(data, 1)
    for counter_bin = 1 : size(data,2)
        length_data(counter_bin_dist, counter_bin) = length(data{counter_bin_dist, counter_bin});
    end
    max_length_data(counter_bin_dist,1) =  max(length_data(counter_bin_dist, :));
    bin_dist{counter_bin_dist,1} = counter_bin_dist *ones( max_length_data(counter_bin_dist,1),1);

    for counter_bin = 1 : size(data,2)
        t_data{counter_bin_dist, counter_bin} = nan(max_length_data(counter_bin_dist,1),1);
        t_data{counter_bin_dist, counter_bin} = data{counter_bin_dist, counter_bin};
        %         t_data{counter_bin_dist, counter_bin}(1:sum(~isnan(data{counter_bin_dist, counter_bin})),1) = data{counter_bin_dist, counter_bin}(~isnan(data{counter_bin_dist, counter_bin}));

    end
end
bin_dist = cell2mat(bin_dist);
t_data = cell2mat(t_data);

% build table
t_data = [bin_dist t_data];
variable_name = cell(1, size(t_data,2)); variable_name{1,1} = 'bin_id';
variable_type =  cell(1, size(t_data,2));  variable_type{1,1} = 'categorical';
for counter_bin = 2 : size(t_data,2)
    variable_name{1,counter_bin} = ['b' num2str(counter_bin-1)];
    variable_type{1,counter_bin} = 'double';
end
table_data = table('Size', [size(t_data,1) size(t_data,2)], 'VariableNames', variable_name, 'VariableTypes', variable_type);

% populate table
for counter_col = 1 : size(table_data,2)
    if counter_col == 1
        table_data(:,counter_col) = table(categorical(t_data(:,counter_col)));
    else
        table_data(:,counter_col) = table(t_data(:,counter_col));
    end
end

bin = (1:(size(table_data,2)-1))';
rm = fitrm(table_data, [variable_name{1,2} '-' variable_name{1,length(variable_name)} ' ~ bin_id'],'WithinDesign',bin);
table_ranova = ranova(rm);
sphericity = mauchly(rm);
corrected = epsilon(rm);
m_ = multcompare(rm,'bin_id');

if (max(double(m_{:,1}) == 3))
    m = [m_{1,3} m_{1,4} m_{1,5}; m_{2,3} m_{2,4} m_{2,5}; m_{4,3} m_{4,4} m_{4,5}];
elseif (max(double(m_{:,1}) == 2))
    m = [m_{1,3} m_{1,4} m_{1,5};];
    % elseif (max(double(m_{:,1}) == 5))
    %     m = [m_{1,3} m_{1,4} m_{1,5}; m_{2,3} m_{2,4} m_{2,5}; m_{4,3} m_{4,4} m_{4,5}];
end
% summary(:,1) = round(m(:,1), 3, 'significant'); summary(:,2) = floor(log10(m(:,2)));
summary = round(m, 3, 'significant');

test = ['F (' num2str(table_ranova{1, 2}), ', ' num2str(table_ranova{3, 2}), ') = ' num2str(round(table_ranova{1, 4}, 2)) ', ' num2str(round(table_ranova{1, 5},3,'significant')) ', ' num2str(round(table_ranova{1, 6},3,'significant')) ', ' num2str(round(table_ranova{1, 7},3,'significant')) newline ...
    'F (' num2str(table_ranova{2, 2}), ', ' num2str(table_ranova{1, 2}), ') = ' num2str(round(table_ranova{2, 4}, 2)) ', ' num2str(round(table_ranova{2, 5},3,'significant')) ', ' num2str(round(table_ranova{2, 6},3,'significant')) ', ' num2str(round(table_ranova{2, 7},3,'significant')) newline ...
    'X^2(' num2str(sphericity{1,3}) ') = ' num2str(round(sphericity{1,2} ,2)) ', ' num2str(round(sphericity{1,4} ,3, 'significant')) ' | ' num2str(round(corrected{1,2} ,3, 'significant'))   ];

% summary = [];

end

%% function build_data_all
function [data_all] = build_data_all(data, bin_dist ,bin)
if isempty(bin_dist)
    flag_time = 1;
else
    flag_time = 0;
end

% determine if work/harvest or sac/lick level
for counter_session = 1 : length(bin)
    length_bin(counter_session) = length(bin{counter_session,1});
end

% build data
if range(length_bin) == 0
    % work/harvest level
    for counter_session = 1 : size(data,1)
        for counter_bin_dist = 1 : max(cell2mat(bin_dist))
            for counter_bin = 1 : max(cell2mat(bin))
                data_ = data{counter_session,1}(1,bin_dist{counter_session,1} == counter_bin_dist &  bin{counter_session,1} == counter_bin );
                data_all{counter_bin_dist, counter_bin}(counter_session,1) = nanmean(data_(~isnan(data_) & data_ ~= 0));
            end
        end
    end
else
    % sac/lick level
    % avg accross time
    if flag_time == 0
        for counter_session = 1 : size(data,1)
            for counter_bin_dist = 1 : max(cell2mat(bin_dist))
                for counter_bin = 1 : max(cell2mat(bin))
                    data_ = data{counter_session,1}(bin_dist{counter_session,1} == counter_bin_dist &  bin{counter_session,1} == counter_bin,: );
                    data_all{counter_bin_dist, counter_bin}(counter_session,1) = nanmean(data_, 'all');
                end
            end
        end
    else
        % preserve time
        for counter_session = 1 : size(data,1)
            for counter_bin = 1 : max(cell2mat(bin))
                for counter_bin_time = 1 : size(data{counter_session,1}, 2)
                    data_ = data{counter_session,1}(bin{counter_session,1} == counter_bin , counter_bin_time);
                    data_all{counter_bin, counter_bin_time}(counter_session,1) = nanmean(data_, 'all');
                end
            end
        end
    end

end
end
%% BUILD FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function extract_population_data
function extract_population_data(path_data_monkey_sorted)
out_path = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep 'session_data' filesep];
mkdir(out_path);

current_monkey_path_ = path_data_monkey_sorted;
fprintf(['### ' 'Analyzing subject ', current_monkey_path_ ' ### \n']);
if ~strcmp(current_monkey_path_(end), filesep)
    current_monkey_path = [current_monkey_path_ filesep];
end
session_list = dir([current_monkey_path filesep '20*' filesep '20*']);
session_list([session_list.isdir] == 0) = [];
sess_list = {session_list.name};

num_sess = length(sess_list);
num_trial_sess = nan(num_sess,1);
num_unit_sess = nan(num_sess,1);
probe_type_sess = nan(num_sess,1);
% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    clearvars EXPERIMENT_DATA LICKS_ALL_DATA SACS_ALL_DATA data_recording
    current_sess = sess_list{counter_sess};
    fprintf(['  ### ' 'sess ', current_sess, ' num. ',...
        num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [current_monkey_path, current_sess(1:7), filesep, current_sess, filesep];

    sess_name = regexprep(current_sess(3:end),'-','');
    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    n_recs = length(rec_list);

    % EXPERIMENT_DATA
    path_to_units = [sess_path 'units' filesep];

    % Load session trial count and electrode type from session meta dat
    session_meta_data_name_ = dir([sess_path '*.xls']);
    session_meta_data_name = session_meta_data_name_.name;
    session_meta_data = readtable([sess_path session_meta_data_name]);

    %EXPERIMENT_DATA.num_trial_sess_rec = (session_meta_data.num_trial);
    EXPERIMENT_DATA.num_trial_sess = sum(session_meta_data.num_trial);
    EXPERIMENT_DATA.probe_type_sess = session_meta_data.elec(1);
    EXPERIMENT_DATA.experiment_type_sess = session_meta_data.exp(1);

    % Count number of units from rec unit summary
    try
        session_rec_unit_summary_name_ = dir([path_to_units '*.mat']);
        session_rec_unit_summary_name = session_rec_unit_summary_name_.name;
        session_rec_unit_summary = load([path_to_units session_rec_unit_summary_name]);
        EXPERIMENT_DATA.num_unit_sess = length(session_rec_unit_summary.cell_ids);
    end
    EXPERIMENT_DATA.sess_date = string(current_sess(3:end));
    EXPERIMENT_DATA.id =  string([num2str(current_sess(3:4)) num2str(current_sess(6:7)) ...
        num2str(current_sess(9:10)) '_combine_' num2str(n_recs)]);
    for counter_rec = 1 : n_recs
        EXPERIMENT_DATA.rec_time(counter_rec, 1) = timeofday(datetime(rec_list{counter_rec}, 'InputFormat', 'yyyy-MM-dd_HH-mm-ss'));

        path_to_rec_eye = [sess_path rec_list{counter_rec} filesep 'analyzed_data' filesep 'behavior_data' filesep 'eye' filesep];
        current_rec_data_eye = dir([path_to_rec_eye '*_ANALYZED_RECAL.mat']);
        current_rec_data_eye_name = current_rec_data_eye.name;
        load([path_to_rec_eye current_rec_data_eye_name], 'TRIALS_DATA');

        EXPERIMENT_DATA.rec_duration(counter_rec, 1) = duration(seconds(TRIALS_DATA.time_end(end) - TRIALS_DATA.time_start(1)));

        path_to_rec_tongue = [sess_path rec_list{counter_rec} filesep 'analyzed_data' filesep 'behavior_data' filesep 'tongue' filesep];
        current_rec_data_tongue_align = dir([path_to_rec_tongue '*_aligned.mat']);
        current_rec_data_tongue_align_name = current_rec_data_tongue_align.name;
        load([path_to_rec_tongue current_rec_data_tongue_align_name], 'align_PD');
        EXPERIMENT_DATA.sample_diff(counter_rec, 1) = align_PD.sample_diff;
        EXPERIMENT_DATA.exp_start_time(counter_rec, 1) = align_PD.BEHAVE_PD_xcorr_time_1K(1,1);

        current_rec_data_tongue = dir([path_to_rec_tongue '*_ANALYZED.mat']);
        current_rec_data_tongue_name = current_rec_data_tongue.name;

        %LICKS_ALL_DATA
        load([path_to_rec_tongue current_rec_data_tongue_name], 'LICKS_ALL_DATA');
        data_recording(counter_rec).LICKS_ALL_DATA.tag = LICKS_ALL_DATA.tag;
        data_recording(counter_rec).LICKS_ALL_DATA.tag_bout = LICKS_ALL_DATA.tag_bout;
        data_recording(counter_rec).LICKS_ALL_DATA.tag_harvest = LICKS_ALL_DATA.tag_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_px = LICKS_ALL_DATA.tongue_tip_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_py = LICKS_ALL_DATA.tongue_tip_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_px = LICKS_ALL_DATA.tongue_mid_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_py = LICKS_ALL_DATA.tongue_mid_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_px = LICKS_ALL_DATA.tongue_r_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_py = LICKS_ALL_DATA.tongue_r_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_px = LICKS_ALL_DATA.tongue_l_px;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_py = LICKS_ALL_DATA.tongue_l_py;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_px_dmax = LICKS_ALL_DATA.tongue_tip_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_tip_py_dmax = LICKS_ALL_DATA.tongue_tip_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_px_dmax = LICKS_ALL_DATA.tongue_mid_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_mid_py_dmax = LICKS_ALL_DATA.tongue_mid_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_px_dmax = LICKS_ALL_DATA.tongue_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_r_py_dmax = LICKS_ALL_DATA.tongue_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_px_dmax = LICKS_ALL_DATA.tongue_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_l_py_dmax = LICKS_ALL_DATA.tongue_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_r_px_dmax = LICKS_ALL_DATA.nose_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_r_py_dmax = LICKS_ALL_DATA.nose_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_l_px_dmax = LICKS_ALL_DATA.nose_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.nose_l_py_dmax = LICKS_ALL_DATA.nose_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_px_onset = LICKS_ALL_DATA.rtube_r_px_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_py_onset = LICKS_ALL_DATA.rtube_r_py_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_px_onset = LICKS_ALL_DATA.rtube_l_px_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_py_onset = LICKS_ALL_DATA.rtube_l_py_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_px_onset = LICKS_ALL_DATA.ltube_r_px_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_py_onset = LICKS_ALL_DATA.ltube_r_py_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_px_onset = LICKS_ALL_DATA.ltube_l_px_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_py_onset = LICKS_ALL_DATA.ltube_l_py_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_px_dmax = LICKS_ALL_DATA.rtube_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_r_py_dmax = LICKS_ALL_DATA.rtube_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_px_dmax = LICKS_ALL_DATA.rtube_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rtube_l_py_dmax = LICKS_ALL_DATA.rtube_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_px_dmax = LICKS_ALL_DATA.ltube_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_r_py_dmax = LICKS_ALL_DATA.ltube_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_px_dmax = LICKS_ALL_DATA.ltube_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.ltube_l_py_dmax = LICKS_ALL_DATA.ltube_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_r_px_dmax = LICKS_ALL_DATA.rew_r_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_r_py_dmax = LICKS_ALL_DATA.rew_r_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_l_px_dmax = LICKS_ALL_DATA.rew_l_px_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_l_py_dmax = LICKS_ALL_DATA.rew_l_py_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_onset = LICKS_ALL_DATA.time_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.time_vmax = LICKS_ALL_DATA.time_vmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_dmax = LICKS_ALL_DATA.time_dmax;
        data_recording(counter_rec).LICKS_ALL_DATA.time_vmin = LICKS_ALL_DATA.time_vmin;
        data_recording(counter_rec).LICKS_ALL_DATA.time_offset = LICKS_ALL_DATA.time_offset;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_dm = LICKS_ALL_DATA.tongue_dm;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm = LICKS_ALL_DATA.tongue_vm;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_ang = LICKS_ALL_DATA.tongue_ang;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_dm_max = LICKS_ALL_DATA.tongue_dm_max;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm_max = LICKS_ALL_DATA.tongue_vm_max;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_vm_min = LICKS_ALL_DATA.tongue_vm_min;
        data_recording(counter_rec).LICKS_ALL_DATA.tongue_ang_max = LICKS_ALL_DATA.tongue_ang_max;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_lick = LICKS_ALL_DATA.duration_lick;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_harvest= LICKS_ALL_DATA.duration_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.duration_harvest= LICKS_ALL_DATA.duration_harvest;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_r_lick_onset= LICKS_ALL_DATA.rew_capacity_r_lick_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_r_lick_offset= LICKS_ALL_DATA.rew_capacity_r_lick_offset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_l_lick_onset= LICKS_ALL_DATA.rew_capacity_l_lick_onset;
        data_recording(counter_rec).LICKS_ALL_DATA.rew_capacity_l_lick_offset= LICKS_ALL_DATA.rew_capacity_l_lick_offset;
        EXPERIMENT_DATA.num_licks_sess_rec(counter_rec, 1) = length(LICKS_ALL_DATA.tag);
        EXPERIMENT_DATA.num_harvests_sess_rec(counter_rec, 1) = sum(LICKS_ALL_DATA.tag_harvest == 1);

        %SACS_ALL_DATA
        load([path_to_rec_eye current_rec_data_eye_name], 'SACS_ALL_DATA');
        data_recording(counter_rec).SACS_ALL_DATA.tag = SACS_ALL_DATA.tag  ;
        data_recording(counter_rec).SACS_ALL_DATA.trial_num = SACS_ALL_DATA.trial_num  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_start = TRIALS_DATA.time_start  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_end = TRIALS_DATA.time_end  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_onset = SACS_ALL_DATA.time_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_vmax = SACS_ALL_DATA.time_vmax  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_offset = SACS_ALL_DATA.time_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.reaction = SACS_ALL_DATA.reaction  ;
        data_recording(counter_rec).SACS_ALL_DATA.duration = SACS_ALL_DATA.duration  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_vm = SACS_ALL_DATA.eye_r_vm  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_amp_m = SACS_ALL_DATA.eye_r_amp_m  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_ang = SACS_ALL_DATA.eye_r_ang  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_px_onset = SACS_ALL_DATA.eye_r_px_onset;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_py_onset = SACS_ALL_DATA.eye_r_py_onset;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_px_onset = SACS_ALL_DATA.visual_px_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_py_onset = SACS_ALL_DATA.visual_py_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_px_offset = SACS_ALL_DATA.eye_r_px_offset;
        data_recording(counter_rec).SACS_ALL_DATA.eye_r_py_offset = SACS_ALL_DATA.eye_r_py_offset;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_px_offset = SACS_ALL_DATA.visual_px_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.tgt_py_offset = SACS_ALL_DATA.visual_py_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.validity = SACS_ALL_DATA.validity  ;
        data_recording(counter_rec).SACS_ALL_DATA.fix_validity = SACS_ALL_DATA.fix_validity  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_fix_offset = SACS_ALL_DATA.time_fix_offset  ;
        data_recording(counter_rec).SACS_ALL_DATA.time_fix_onset = SACS_ALL_DATA.time_fix_onset  ;
        data_recording(counter_rec).SACS_ALL_DATA.x_fix_after = SACS_ALL_DATA.x_fix_after  ;
        data_recording(counter_rec).SACS_ALL_DATA.y_fix_after = SACS_ALL_DATA.y_fix_after  ;
        data_recording(counter_rec).SACS_ALL_DATA.x_fix_before = SACS_ALL_DATA.x_fix_before  ;
        data_recording(counter_rec).SACS_ALL_DATA.y_fix_before = SACS_ALL_DATA.y_fix_before;
        for counter_trial = 1:length(TRIALS_DATA.time_state_cue_present)
            data_recording(counter_rec).SACS_ALL_DATA.num_trial_attempt(counter_trial) = length(TRIALS_DATA.time_state_cue_present{counter_trial})  ;
        end

        EXPERIMENT_DATA.num_sacs_sess_rec(counter_rec, 1) = length(SACS_ALL_DATA.tag);
        EXPERIMENT_DATA.num_trial_sess_rec(counter_rec, 1) = length(TRIALS_DATA.time_start);
        EXPERIMENT_DATA.num_fixs_sess_rec(counter_rec, 1) = length(SACS_ALL_DATA.time_fix_onset);

        %PUPILS_DATA
        load([path_to_rec_eye current_rec_data_eye_name], 'PUPILS_DATA');
        data_recording(counter_rec).PUPILS_ALL_DATA.pupil_area = PUPILS_DATA.eye_r_pupil_a_filt;
        data_recording(counter_rec).PUPILS_ALL_DATA.time_1K = PUPILS_DATA.time_1K;
        data_recording(counter_rec).PUPILS_ALL_DATA.validity = PUPILS_DATA.validity;
        EXPERIMENT_DATA.length_pupil_sess_rec(counter_rec, 1) = length(PUPILS_DATA.eye_r_pupil_a_filt);

        %BLINKS_DATA
        load([path_to_rec_eye current_rec_data_eye_name], 'BLINKS_DATA');
        data_recording(counter_rec).BLINKS_ALL_DATA.blink_onset = BLINKS_DATA.eye_r_blink_onset;
        data_recording(counter_rec).BLINKS_ALL_DATA.blink_offset = BLINKS_DATA.eye_r_blink_offset;
        EXPERIMENT_DATA.num_blink_sess_rec(counter_rec, 1) = length(BLINKS_DATA.eye_r_blink_onset);

        path_to_raw = [path_to_rec_tongue '..' filesep '..' filesep '..' filesep 'raw_data' filesep];

        % add vid frame to EXPERIMENT_DATA
        % Save data
        if ~isempty(dir([path_to_rec_tongue '*_DLC.mp4']))
            vid = dir([path_to_rec_tongue '*_DLC.mp4']);
        elseif ~isempty(dir([path_to_raw '*.mp4']))
            vid = dir([path_to_raw '*.mp4']);
        else
            vid = dir([path_to_raw '*.avi']);
        end

        vid_name = vid.name;
        vid_path = vid.folder;

        EXPERIMENT_DATA.vid_frame = readFrame(VideoReader([vid_path filesep vid_name]));
    end

    % compute break duration -> time between recordings
    for counter_rec = 1 : n_recs - 1
        EXPERIMENT_DATA.break_duration(counter_rec,1) = EXPERIMENT_DATA.rec_time(counter_rec + 1) - ...
            EXPERIMENT_DATA.rec_time(counter_rec) - EXPERIMENT_DATA.rec_duration(counter_rec);
    end
    EXPERIMENT_DATA.break_duration(n_recs,1) = 0;

    % concatenate LICKS_ALL_DATA, SACS_ALL_DATA, BLINKS_ALL_DATA, PUPILS_ALL_DATA
    [LICKS_ALL_DATA,SACS_ALL_DATA,BLINKS_ALL_DATA,PUPILS_ALL_DATA] = PGH_combine_recs(data_recording);

    % add num_licks and num_sacs to EXPERIMENT_DATA
    EXPERIMENT_DATA.num_licks_sess = length(LICKS_ALL_DATA.tag);
    EXPERIMENT_DATA.num_sacs_sess = length(SACS_ALL_DATA.tag);

    save([out_path  char(EXPERIMENT_DATA.id) '.mat'], 'LICKS_ALL_DATA','SACS_ALL_DATA','BLINKS_ALL_DATA','PUPILS_ALL_DATA','EXPERIMENT_DATA','-v7.3')

end
fprintf('### ALL DONE. ###\n')
end
%% function build_population_data
function build_population_data(path_data_monkey_sorted)
%% Load data
path_BEHAVE_session_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep 'session_data'];
out_path = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted filesep];

if ~strcmp(path_BEHAVE_session_data(end), filesep);path_BEHAVE_session_data = [path_BEHAVE_session_data filesep];end
session_list = dir([path_BEHAVE_session_data '*mat']);
num_session = numel(session_list);

%% Loop over sessions
for counter_session = 1 : num_session
    fprintf(['### ' 'Building population data session num. ' num2str(counter_session), ' / ' num2str(num_session) ' ### \n']);
    load([session_list(counter_session).folder filesep session_list(counter_session).name]);

    % population_experiment_data
    fieldnames_EXPERIMENT_DATA = fieldnames(EXPERIMENT_DATA);
    for counter_fieldnames = 1 : length(fieldnames_EXPERIMENT_DATA)
        population_experiment_data.(fieldnames_EXPERIMENT_DATA{counter_fieldnames}){counter_session,1} = EXPERIMENT_DATA.(fieldnames_EXPERIMENT_DATA{counter_fieldnames});
    end

    % population_lick_data
    fieldnames_LICK_DATA = fieldnames(LICKS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_LICK_DATA)
        population_lick_data.(fieldnames_LICK_DATA{counter_fieldnames}){counter_session,1} = LICKS_ALL_DATA.(fieldnames_LICK_DATA{counter_fieldnames});
    end

    % population_sac_data
    fieldnames_SAC_DATA = fieldnames(SACS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_SAC_DATA)
        population_sac_data.(fieldnames_SAC_DATA{counter_fieldnames}){counter_session,1} = SACS_ALL_DATA.(fieldnames_SAC_DATA{counter_fieldnames});
    end

    % population_blink_data
    fieldnames_BLINK_DATA = fieldnames(BLINKS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_BLINK_DATA)
        population_blink_data.(fieldnames_BLINK_DATA{counter_fieldnames}){counter_session,1} = BLINKS_ALL_DATA.(fieldnames_BLINK_DATA{counter_fieldnames});
    end

    % population_pupil_data
    fieldnames_PUPIL_DATA = fieldnames(PUPILS_ALL_DATA);
    for counter_fieldnames = 1 : length(fieldnames_PUPIL_DATA)
        population_pupil_data.(fieldnames_PUPIL_DATA{counter_fieldnames}){counter_session,1} = PUPILS_ALL_DATA.(fieldnames_PUPIL_DATA{counter_fieldnames});
    end
end
%% Save data
fprintf(['Saving .mat file' ' ...'])
save([out_path 'population_experiment_data' '.mat'], 'population_experiment_data', '-v7.3');
save([out_path 'population_lick_data' '.mat'], 'population_lick_data', '-v7.3');
save([out_path 'population_sac_data' '.mat'], 'population_sac_data', '-v7.3');
save([out_path 'population_blink_data' '.mat'], 'population_blink_data', '-v7.3');
save([out_path 'population_pupil_data' '.mat'], 'population_pupil_data', '-v7.3');

%rmdir(path_BEHAVE_session_data)
fprintf('### ALL DONE. ###\n')

end
%% function generate_population_meta_data
function generate_population_meta_data(path_data_monkey_sorted)
%% Load data
path_BEHAVE_population_data = ['D:' filesep 'BEHAVE' filesep path_data_monkey_sorted ];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
out_path = [path_BEHAVE_population_data];
load([path_BEHAVE_population_data 'population_experiment_data.mat']);
%% Build inclusion table
num_session = length(population_experiment_data.sess_date);
include_tongue = ones(num_session,1);
include_eye = ones(num_session,1);

session_date = population_experiment_data.sess_date;
num_trial = population_experiment_data.num_trial_sess;
num_lick = population_experiment_data.num_licks_sess;
num_sac = population_experiment_data.num_sacs_sess;
exp_type = population_experiment_data.experiment_type_sess;

meta_data = table(session_date, exp_type, num_trial, num_lick, num_sac, include_tongue, include_eye);

%% Save
writetable(meta_data, [out_path 'meta_data.xls'])
end
%% function build_foraging_analysis
function build_foraging_analysis(params,path_data_monkey_sorted)
%% Load data
fprintf('Loading data ... ')
path.path_data_monkey_sorted = path_data_monkey_sorted;
path_BEHAVE_population_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
path.out_path = [path_BEHAVE_population_data];
load([path_BEHAVE_population_data 'population_experiment_data.mat']); % exp data
load([path_BEHAVE_population_data 'population_lick_data.mat']); % lick data
load([path_BEHAVE_population_data 'population_sac_data.mat']); % sac data
load([path_BEHAVE_population_data 'population_blink_data.mat']); % blink data
load([path_BEHAVE_population_data 'population_pupil_data.mat']); % pupil data
meta_data = readtable([path_BEHAVE_population_data 'meta_data.xls']); % meta data
include_tongue = find(meta_data{:,6} == 1);
include_eye = find(meta_data{:,7} == 1);
%include = find(meta_data{:,13} == 1);

% include = include_tongue(ismember(include_tongue,include_eye));
include = 1:size(meta_data,1);

% add weight data
for counter_session = 1 : length(population_experiment_data.sess_date)
    population_experiment_data.weight(counter_session,1) = table2cell(meta_data(counter_session,8));
end

% exp data
field_name = string(fieldnames(population_experiment_data));
for counter_field_name = 1 : length(field_name)
    population_experiment_data.(field_name(counter_field_name)) = population_experiment_data.(field_name(counter_field_name))(include,1);
end

% lick data
field_name = string(fieldnames(population_lick_data));
for counter_field_name = 1 : length(field_name)
    population_lick_data.(field_name(counter_field_name)) = population_lick_data.(field_name(counter_field_name))(include,1);
end

% sac data
field_name = string(fieldnames(population_sac_data));
for counter_field_name = 1 : length(field_name)
    population_sac_data.(field_name(counter_field_name)) = population_sac_data.(field_name(counter_field_name))(include,1);
end

% blink data
field_name = string(fieldnames(population_blink_data));
for counter_field_name = 1 : length(field_name)
    population_blink_data.(field_name(counter_field_name)) = population_blink_data.(field_name(counter_field_name))(include,1);
end

% pupil data
field_name = string(fieldnames(population_pupil_data));
for counter_field_name = 1 : length(field_name)
    population_pupil_data.(field_name(counter_field_name)) = population_pupil_data.(field_name(counter_field_name))(include,1);
end
fprintf('---> Complete. \n')

%% Build variables
fprintf('Building variables ... ')
[count, lick, sac] = build_variables(population_experiment_data, population_lick_data, population_sac_data, population_blink_data, population_pupil_data, path);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% compute model
fprintf('Building model ... ')
[lick] = build_model(lick,sac,count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')
%% compute endpoint error mag
fprintf('Computing endpoint error magnitude ... ')
[lick] = endpoint_error_mag(lick,count);
[sac] = endpoint_error_mag(sac,count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% pad/cut variables
fprintf('Padding data ... ')
[lick] = pad_data(lick, count); % pad lick data
[sac] = pad_data(sac, count); % pad sac data

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% build bins
fprintf('Building bins ... ')
[count] = build_bins(lick, sac, count, path);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')
%% add RPE
fprintf('Add RPE  ... ')
[count] = add_RPE(lick, count);

clearvars -except params lick sac path count
fprintf('---> Complete. \n')

%% compute move type prob
% fprintf('Computing move type probability... ')
% [lick] = move_type_prob(lick,count);
% [sac] = move_type_prob(sac,count);
%
% clearvars -except params lick sac path count
% fprintf('---> Complete. \n')

%% save data structures
fprintf('Saving... ')
save([path.out_path  'lick.mat' ], 'lick', '-v7.3');
save([path.out_path  'sac.mat' ], 'sac', '-v7.3');
save([path.out_path  'count.mat' ], 'count', '-v7.3');
clearvars -except params lick sac path count
fprintf('---> Complete. \n')
end
%% PLOT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function plot_session_analysis
function plot_session_analysis(params,path_data_monkey_sorted)
%% Load data
path_BEHAVE_population_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
out_path = [path_BEHAVE_population_data filesep 'SESSION' filesep];
load([path_BEHAVE_population_data 'population_experiment_data.mat']);
load([path_BEHAVE_population_data 'count']);
meta_data = readtable([path_BEHAVE_population_data 'meta_data.xls']); % meta data
include_tongue = find(meta_data{:,6} == 1);
include_eye = find(meta_data{:,7} == 1);

include = include_tongue(ismember(include_tongue,include_eye));
%% Build data
probe_type_sess = cell2mat(population_experiment_data.probe_type_sess);
num_trial_sess = cell2mat(population_experiment_data.num_trial_sess);
num_unit_sess   = cell2mat(population_experiment_data.num_unit_sess);
num_licks = cell2mat(population_experiment_data.num_licks_sess);
num_sacs = cell2mat(population_experiment_data.num_sacs_sess);

mean_num_sacs = nanmean(num_sacs);
std_num_sacs = nanstd(num_sacs);

mean_num_licks = nanmean(num_licks);
std_num_licks = nanstd(num_licks);

mean_num_unit_sess = nanmean(num_unit_sess);
std_num_unit_sess = nanstd(num_unit_sess);

% silicon array
is_sa = probe_type_sess == 1 | probe_type_sess == 2;
num_sess_sa = sum(is_sa);

num_trial_sa = sum(num_trial_sess(is_sa));
mean_num_trial_sa = nanmean(num_trial_sess(is_sa));
std_num_trial_sa = nanstd(num_trial_sess(is_sa));
se_num_trial_sa = std_num_trial_sa/sqrt(num_sess_sa);

try
    num_unit_sa = sum(num_unit_sess(is_sa));
    mean_num_unit_sa = nanmean(num_unit_sess(is_sa));
    std_num_unit_sa = nanstd(num_unit_sess(is_sa));
    se_num_unit_sa = std_num_unit_sa/sqrt(num_sess_sa);
catch
    num_unit_sa = 0;
    mean_num_unit_sa = 0;
    std_num_unit_sa = 0;
    se_num_unit_sa = 0;
end

if num_unit_sa == 0
    mean_num_unit_sa = 0;
    std_num_unit_sa = 0;
    se_num_unit_sa = 0;
end

% hept/tet
is_ht = probe_type_sess == 3 | probe_type_sess == 4;
num_sess_ht = sum(is_ht);

num_trial_ht = sum(num_trial_sess(is_ht));
mean_num_trial_ht = nanmean(num_trial_sess(is_ht));
std_num_trial_ht = nanstd(num_trial_sess(is_ht));
se_num_trial_ht = std_num_trial_ht/sqrt(num_sess_ht);

try
    num_unit_ht = sum(num_unit_sess(is_ht));
    mean_num_unit_ht = nanmean(num_unit_sess(is_ht));
    std_num_unit_ht = nanstd(num_unit_sess(is_ht));
    se_num_unit_ht = std_num_unit_ht/sqrt(num_sess_ht);
catch
    num_unit_ht = 0;
    mean_num_unit_ht = 0;
    std_num_unit_ht = 0;
    se_num_unit_ht = 0;
end
% combinded
num_sess_all = num_sess_ht + num_sess_sa;

num_trial_all = num_trial_ht + num_trial_sa;
mean_num_trial_all = nanmean(num_trial_sess);
std_num_trial_all = nanstd(num_trial_sess);
se_num_trial_all = std_num_trial_all/sqrt(num_sess_all);

num_unit_all = num_unit_ht + num_unit_sa;

% extract num_unit_sess_type
num_unit_sess_sa = num_unit_sess;
num_unit_sess_sa(~is_sa) = nan;

num_unit_sess_ht = num_unit_sess;
num_unit_sess_ht(~is_ht) = nan;

for counter_session = 1 : num_sess_all
    rec_duration(counter_session, 1) = sum(population_experiment_data.rec_duration{counter_session, 1});
    break_duration(counter_session, 1) = sum(population_experiment_data.break_duration{counter_session, 1});
    session_duration(counter_session, 1) = rec_duration(counter_session, 1) + break_duration(counter_session, 1);
end

mean_rec_duration = mean(rec_duration);
std_rec_duration = std(rec_duration);

mean_break_duration = mean(break_duration);
std_break_duration = std(break_duration);

mean_session_duration = mean(session_duration);
std_session_duration = std(session_duration);

rec_duration.Format = 'hh:mm:ss';
break_duration.Format = 'hh:mm:ss';
session_duration.Format = 'hh:mm:ss';

%% Plot data - Fig 1
x_axis_session = 1:num_sess_all;
figure
subplot(3,3,1)
bar(x_axis_session, session_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('session duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,2)
bar(x_axis_session, rec_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('rec duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,3)
bar(x_axis_session, break_duration, 'FaceColor', 'black')
xlabel('recording session')
ylabel('time')
title('break duration')
xlim([0 length(num_trial_sess)+1])
ylim([duration([0 0 0]) duration([4 0 0])])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,4)
bar(x_axis_session, num_trial_sess, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. correct trials')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,5)
bar(x_axis_session, num_sacs, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. sacs')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,6)
bar(x_axis_session, num_licks, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. licks')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,7)
%bar(x_axis_session, num_unit_sess_sa, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. neurons: silicon array')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
xlim([0 length(num_trial_sess)+1])
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,8)
%bar(x_axis_session, num_unit_sess_ht, 'FaceColor', 'black')
xlabel('recording session')
ylabel('count')
title('num. neurons: tetrode/heptode')
xlim([0 length(num_trial_sess)+1])
xticks(1:1:length(num_trial_sess))
set(gca,'xticklabel',population_experiment_data.sess_date  ,'fontsize',1, 'TickLabelInterpreter','none')

subplot(3,3,9)
pie([num_sess_sa, num_sess_ht], {['silicon array: ' num2str(num_sess_sa)], ['terode/heptode: ' num2str(num_sess_ht)]})
title('num. sessions w/ probe type')

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [out_path 'session_list_summary'], 'pdf')

%% Plot data - Fig 2
figure
subplot(3,3,1)
histogram(session_duration, 'FaceColor', 'k');
xline(mean_session_duration, 'r', 'LineWidth', 1)
xline(mean_session_duration + std_session_duration, '--r', 'LineWidth', 0.5)
xline(mean_session_duration - std_session_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('sess. duration')

subplot(3,3,2)
histogram(rec_duration, 'FaceColor', 'k');
xline(mean_rec_duration, 'r', 'LineWidth', 1)
xline(mean_rec_duration + std_rec_duration, '--r', 'LineWidth', 0.5)
xline(mean_rec_duration - std_rec_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('rec. duration')

subplot(3,3,3)
histogram(break_duration, 'FaceColor', 'k');
xline(mean_break_duration, 'r', 'LineWidth', 1)
xline(mean_break_duration + std_break_duration, '--r', 'LineWidth', 0.5)
xline(mean_break_duration - std_break_duration, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('break duration')

subplot(3,3,4)
histogram(num_trial_sess, 'FaceColor', 'k');
xline(mean_num_trial_all, 'r', 'LineWidth', 1)
xline(mean_num_trial_all + std_num_trial_all, '--r', 'LineWidth', 0.5)
xline(mean_num_trial_all - std_num_trial_all, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. correct trials')

subplot(3,3,5)
histogram(num_sacs, 'FaceColor', 'k');
xline(mean_num_sacs, 'r', 'LineWidth', 1)
xline(mean_num_sacs + std_num_sacs, '--r', 'LineWidth', 0.5)
xline(mean_num_sacs - std_num_sacs, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. sacs')

subplot(3,3,6)
histogram(num_licks, 'FaceColor', 'k');
xline(mean_num_licks, 'r', 'LineWidth', 1)
xline(mean_num_licks + std_num_licks, '--r', 'LineWidth', 0.5)
xline(mean_num_licks - std_num_licks, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. licks')

subplot(3,3,7)
%histogram(num_unit_sess(is_ht), 'FaceColor', 'k');
xline(mean_num_unit_ht, 'r', 'LineWidth', 1)
xline(mean_num_unit_ht + std_num_unit_ht, '--r', 'LineWidth', 0.5)
xline(mean_num_unit_ht - std_num_unit_ht, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. units w/ hept/tet')

subplot(3,3,8)
%histogram(num_unit_sess(is_sa), 'FaceColor', 'k');
xline(mean_num_unit_sa, 'r', 'LineWidth', 1)
xline(mean_num_unit_sa + std_num_unit_sa, '--r', 'LineWidth', 0.5)
xline(mean_num_unit_sa - std_num_unit_sa, '--r', 'LineWidth', 0.5)
ylabel('count')
xlabel('num. units w/ silicon array')

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [out_path 'session_hist_summary'], 'pdf')
%% Plot data - Fig 3
figure
subplot(3,3,1)
plot(rec_duration,num_trial_sess, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. correct trials')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,2)
plot(rec_duration,num_sacs, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. sacs')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,3)
plot(rec_duration,num_licks, '.k', 'MarkerSize', 20);
xlabel('rec. duration')
ylabel('num. licks')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,4)
plot(break_duration,num_trial_sess, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. correct trials')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,5)
plot(break_duration,num_sacs, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. sacs')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,6)
plot(break_duration,num_licks, '.k', 'MarkerSize', 20);
xlabel('break duration')
ylabel('num. licks')
ylim([0 inf])
xlim([duration([0 0 0]) duration([3 0 0])])

subplot(3,3,7)
plot(num_trial_sess,num_sacs, '.k', 'MarkerSize', 20);
xlabel('num. correct trials')
ylabel('num. sacs')
ylim([0 inf])
xlim([0 inf])

subplot(3,3,8)
plot(num_trial_sess,num_licks, '.k', 'MarkerSize', 20);
xlabel('num. correct trials')
ylabel('num. licks')
ylim([0 inf])
xlim([0 inf])

subplot(3,3,9)
plot(num_sacs,num_licks, '.k', 'MarkerSize', 20);
xlabel('num. sacs')
ylabel('num. licks')
ylim([0 inf])
xlim([0 inf])

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [out_path 'session_list_comparison'], 'pdf')
%% Plot data - Fig 4
figure;
clear session_list
session_list = population_experiment_data.sess_date;

num_session = length(session_list);

num_col = floor(num_session/5);
num_row = ceil(num_session/floor(num_session/5));

for counter_session = 1 : num_session
    included = [];
    if ismember(counter_session, include)
        included = '*';
    end
    subplot(num_row, num_col, counter_session)
    imshow(population_experiment_data.vid_frame{counter_session, 1}  )
    title([char(session_list{counter_session,1}) included], 'Interpreter', 'none')
end

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all)  ], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [out_path 'session_video_check'], 'pdf')

%% Plot data - Fig 5
figure;
clear session_list
for counter_session = 1 : length(population_experiment_data.sess_date)
    session_list(counter_session) = (population_experiment_data.sess_date{counter_session});
end

session_list = session_list(include);
for counter_session = 1 : count.num_session
    session_list(counter_session) = (population_experiment_data.sess_date{counter_session});
    weight = count.weight_bin{counter_session,1}(1,1);
    dist_tube_bin_l = mode(count.dist_tube_bin{counter_session,1}(count.dir_harvest_bin{counter_session, 1}  == 2));
    dist_tube_bin_r = mode(count.dist_tube_bin{counter_session,1}(count.dir_harvest_bin{counter_session, 1}  == 1));

    subplot(num_row, num_col, counter_session)
    x_axis = categorical({'L', 'R'}); x_axis = reordercats(x_axis, {'L', 'R'});
    bar(x_axis, [count.dist_l_tube_harvest_mean(counter_session,1), count.dist_r_tube_harvest_mean(counter_session,1)],'FaceColor', 'k')
    title([char(session_list(counter_session)) ' | l:' num2str(dist_tube_bin_l) ', r:' num2str(dist_tube_bin_r)  ', w: ' num2str(weight)])
    xlabel('tube')
    ylabel('tube dist (mm)')
    ylim([0 13])
end

sgtitle([path_data_monkey_sorted ' | sessions: ' num2str(num_sess_all) ], 'Interpreter', 'none')
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [out_path 'session_tube_dist_weight' num2str(get(gcf).Number)], 'pdf')



end
%% function plot_foraging_analysis
function plot_foraging_analysis(params,path_data_monkey_sorted)
%% load data structures
fprintf('Loading... ')
path.path_data_monkey_sorted = path_data_monkey_sorted;
path_BEHAVE_population_data = [pwd  filesep 'BEHAVE' filesep path_data_monkey_sorted];
if ~strcmp(path_BEHAVE_population_data(end), filesep);path_BEHAVE_population_data = [path_BEHAVE_population_data filesep];end
path.out_path = [path_BEHAVE_population_data];
lick = []; sac = []; count = [];

load([path.out_path  'lick.mat'],'lick' );
load([path.out_path  'sac.mat' ], 'sac');
load([path.out_path  'count.mat' ],'count');
colors = ['k'; 'r'; 'g'; 'b'; 'c'; 'm';'y'; 'k'; 'r'; 'g'; 'b'; 'c'; 'm';'y'];

meta_data = readtable([path_BEHAVE_population_data 'meta_data.xls']); % meta data

count.disclude_work_r = find(meta_data{:,9} ~= 1);
count.disclude_sac_r = find(meta_data{:,10} ~= 1);
count.disclude_harvest_r = find(meta_data{:,11} ~= 1);
count.disclude_lick_r = find(meta_data{:,12} ~= 1);
count.disclude_work_l = find(meta_data{:,13} ~= 1);
count.disclude_sac_l = find(meta_data{:,14} ~= 1);
count.disclude_harvest_l = find(meta_data{:,15} ~= 1);
count.disclude_lick_l = find(meta_data{:,16} ~= 1);

plot_session = 0;
if plot_session == 1
    count_session = count.num_session;
else
    count_session = 1;
    ext = [];
end
fprintf('---> Complete. \n')

%% Figure - work period
%%% PLOT (1)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
figure
hold on
for counter_session = 1 : count.num_session
    mode_data(counter_session,1) = nanmean(count.dist_tube_bin{counter_session, 1});
end
x_axis = categorical({'near', 'mid', 'far'}); x_axis = reordercats(x_axis,{'near', 'mid', 'far'});
data = count.num_trial;
for counter_bin = 1 : max(mode_data)
    data_{counter_bin} = data(mode_data==counter_bin);
    mean_data(counter_bin) = nanmean(data_{counter_bin});
    sem_data(counter_bin) = nanstd(data_{counter_bin})/sqrt(length(data_{counter_bin}));
    %     plot(x_axis(counter_bin),data_{counter_bin},'or','MarkerSize', 5)
end
errorbar(x_axis, mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
plot(x_axis, mean_data, colors(1), 'LineWidth',1)
xlabel('tube dist');
ylabel('# completed trials in session')
title([num2str(mean_data(1)) '+/-' num2str(sem_data(1)) ' | ' num2str(mean_data(2)) '+/-' num2str(sem_data(2)) ' | ' num2str(mean_data(3)) '+/-' num2str(sem_data(3))])
ESN_Beautify_Plot(gcf, [5,5])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_trials' ], 'pdf');
%close(gcf)

%%% PLOT (2)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["num_trial_work", "num_trial_attempt_work" "duration_work", "rew_gained_work"];
figure
for counter_variable = 1 : length(variable_list)
    data = cell2mat(sac.(variable_list(counter_variable)));data = reshape(data', numel(data),1);
    % for osf           
    writetable(array2table(data), [path.out_path 'data' filesep path.path_data_monkey_sorted '_work_hist_' char(variable_list(counter_variable))  ext '.csv']);

    mean_data = nanmean(data) ;
    std_data = nanstd(data)/sqrt(sum(~isnan(data)));
    if contains(variable_list(counter_variable), 'num_trial_work')
        x_lim = [0 20];
        x_ticks = [0:5:20];
        data(data<1 | data > 20) = nan;
    elseif contains(variable_list(counter_variable), 'num_trial_attempt_work')
        x_lim = [0 40];
        x_ticks = [0:5:40];
        data(data<1 | data > 40) = nan;
    elseif contains(variable_list(counter_variable), 'duration_work')
        x_lim = [0 60];
        x_ticks = [0:10:60];
        data(data<1 | data > 60) = nan;
    elseif contains(variable_list(counter_variable), 'rew_gained_work')
        x_lim = [0 1];
        x_ticks = [0:0.2:1];
    end
    subplot(1,length(variable_list),counter_variable)
    hold on
    histogram(data,'FaceColor', 'k','Normalization','probability')
    xline(mean_data,'-r', 'LineWidth',2)
    xline(mean_data+std_data,'--r', 'LineWidth',1)
    xline(mean_data-std_data,'--r', 'LineWidth',1)
    xlabel(variable_list(counter_variable),'Interpreter', 'none')
    ylabel('prob.')
    xlim(x_lim)
    xticks(x_ticks)
    title([num2str(mean_data) '+/-' num2str(std_data)])
end
sgtitle(['Total trials: ' num2str(nanmean(count.num_trial)) ' +/- ' num2str(nanstd(count.num_trial)/sqrt(length(count.num_trial))) ' | Total intervals: ' num2str(nanmean(count.num_harvest))  ' +/- ' num2str(nanstd(count.num_harvest)/sqrt(length(count.num_harvest)))])
ESN_Beautify_Plot(gcf, [12,3])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_work_hist' ], 'pdf');
%close(gcf)

%%% PLOT (3)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["num_trial_work", "num_trial_attempt_work", "duration_work", "rew_gained_work"];
bin_list = ["interval_bin", "dist_tube_bin", "weight_bin"];
for counter_session = 1 : count_session
    figure
    for counter_bin_list = 1 : length(bin_list)
        for counter_variable = 1 : length(variable_list)
            mean_data = [];
            sem_data = [];
            if counter_bin_list==1
                shift = 0;
            else
                shift = (counter_bin_list - 1)*length(variable_list);
            end
            y_lim = ([-inf inf]);
            y_ticks = 'auto';
            subplot(length(bin_list),length(variable_list),counter_variable + shift)
            hold on
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_work_l, count.disclude_work_r, count.dir_harvest_bin);
            end
            bin = cell2mat(bin);
            bin_dir = cell2mat(count.dir_harvest_bin); is_bin_dir = bin_dir == 1 | bin_dir == 2;
            for counter_bin = 1 : max(bin)
                data = cell2mat(sac.(variable_list(counter_variable))); data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data) | data<=0)=[];
                data_all{1,counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
            end
            % compute ANOVA
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_work_Xbin' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            % build
            data_all = [];
            errorbar(1:max(bin), mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(1), 'LineWidth',1)
            xlabel(bin_list(counter_bin_list),'Interpreter','none');
            ylabel(variable_list(counter_variable),'Interpreter', 'none')
            xlim([0.5 max(bin)+0.5])
            xticks([1:max(bin)])
            ylim(y_lim)
            yticks(y_ticks)
            title(test, 'Interpreter', 'none')
            text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
            if plot_session == 1
                % to plot session
                mean_data = [];
                sem_data = [];
                bin = count.(bin_list(counter_bin_list)){counter_session,1};
                bin_dir = count.dir_harvest_bin{counter_session,1}; is_bin_dir = bin_dir == 1 | bin_dir == 2;
                for counter_bin = 1 : max(bin)
                    data = sac.(variable_list(counter_variable)){counter_session,1}; data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data) | data<=0)=[];
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                    %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
                end
                errorbar(1:max(bin), mean_data, sem_data, colors(2), 'LineWidth', 2, 'MarkerSize', 10);
                %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
                plot(1:max(bin), mean_data, colors(2), 'LineWidth',1)
                ext = [ '_' num2str(counter_session)];
                path_ext =  ['SESSION' filesep 'work' filesep ];
            end
        end
    end
    ESN_Beautify_Plot(gcf, [20,10])
    saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_work_Xbin' ext], 'pdf');
    %close(gcf)
end
%%% PLOT (3.1) %%%
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_work_l, count.disclude_work_r, count.dir_harvest_bin_);
        end
        [bin_dist] = disclude_data(count.dist_tube_bin, count.disclude_work_l, count.disclude_work_r, count.dir_harvest_bin);

        %         % for RANOVA
        %         [data_all] = build_data_all(sac.(variable_list(counter_variable)), bin_dist ,bin);
        %         % compute ANOVA
        %         if ~contains(bin_list(counter_bin_list),'dist')
        %             [test, summary] = PGH_stat_RANOVA(data_all);
        %             data_all = [];
        %         else
        %             test = [];
        %             summary = [];
        %         end

        bin = cell2mat(bin);
        bin_dist = cell2mat(bin_dist);
        bin_dir = cell2mat(count.dir_harvest_bin); is_bin_dir = bin_dir == 1 | bin_dir == 2;

        for counter_bin_dist = 1 : max(bin_dist)
            for counter_bin = 1 : max(bin)
                data = cell2mat(sac.(variable_list(counter_variable))); data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:); data(isnan(data) | data<=0)=[];
                data_all{counter_bin_dist, counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_bin_dist), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(counter_bin_dist), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_bin_dist), 'LineWidth',1)
        end
        % compute ANOVA
        [test, summary, table_data] = PGH_stat_ANOVA(data_all);
        % for osf
        writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_work_Xbin_xdist' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

        data_all = [];
        xlabel(bin_list(counter_bin_list),'Interpreter','none');
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
        title(test, 'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_work_Xbin_xdist' ext], 'pdf');
%close(gcf)

%%% PLOT (4)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list_x = ["pupil_area","pupil_area", "eye_vigor_w",  "eye_vigor_w"];
variable_list_y = ["eye_vigor_w","duration_work",  "dist_err_eye_tgt","var_err_eye_tgt"];
bin_list = ["combined", "interval_bin", "dist_tube_bin", "weight_bin"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list_x)
        mean_data_y_bin = [];
        sem_data_y_bin = [];
        if contains(variable_list_x(counter_variable), 'num_trial_work') %&& contains(bin_list(counter_bin_list), 'combined')
            bins = 1 : 10;
            edges = 0.5 : 1 : 10.5;
            x_ticks = [bins];
        elseif contains(variable_list_x(counter_variable), 'pupil_area') %&& contains(bin_list(counter_bin_list), 'combined')
            %             if contains(path.path_data_monkey_sorted,'59d')
            %                 bins = 5500 :250: 7750;
            %                 edges = 5250 : 250 : 8000;
            %                 x_ticks = [bins];
            %             elseif contains(path.path_data_monkey_sorted,'125d')
            %                 bins = 5000 :250: 7000;
            %                 edges = 4750 : 250 : 7250;
            %                 x_ticks = [bins];
            %             end
            bins = -1 :0.2: 1;
            edges = -1.1 : 0.2 : 1.1;
            x_ticks = [bins];
        elseif contains(variable_list_x(counter_variable), 'eye_vigor_w') %&& contains(bin_list(counter_bin_list), 'combined')
            bins = 0.6 :0.1: 1.4;
            edges = 0.5 : 0.1 : 1.5;
            x_ticks = [bins];
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list_x);
        end
        subplot(length(bin_list) ,length(variable_list_x),counter_variable + shift)
        data_x = cell2mat(sac.(variable_list_x(counter_variable)));
        if ~contains(variable_list_y(counter_variable),'var_err_eye_tgt')
            data_y = cell2mat(sac.(variable_list_y(counter_variable)));
        else
            data_y = cell2mat(sac.err_eye_tgt_py_rot);
            data_x_y = cell2mat(sac.err_eye_tgt_px_rot);
            data_y_y = cell2mat(sac.err_eye_tgt_py_rot);
        end
        if size(data_x,1) == count.num_session && size(data_y,1) == count.num_session % both harvest level
            data_x = reshape(data_x', numel(data_x),1);
            if ~contains(variable_list_y(counter_variable),'var_err_eye_tgt')
                data_y = reshape(data_y', numel(data_y),1);
            else
                data_x_y = reshape(data_x_y', numel(data_x_y),1);
                data_y_y = reshape(data_y_y', numel(data_y_y),1);
            end
            flag_harvest_level_x = 1;
            flag_harvest_level_y = 1;
        elseif size(data_x,1) ~= count.num_session && size(data_y,1) ~= count.num_session % both lick level
            flag_harvest_level_x = 0;
            flag_harvest_level_y = 0;
        elseif size(data_x,1) ~= count.num_session && size(data_y,1) == count.num_session % x lick level, y harvest  level
            data_x = nanmean(data_x,2);
            if ~contains(variable_list_y(counter_variable),'var_err_eye_tgt')
                data_y = reshape(data_y', numel(data_y),1);data_y(isnan(data_y)) = [];
            else
                data_x_y = reshape(data_x_y', numel(data_x_y),1);data_x_y(isnan(data_x_y)) = [];
                data_y_y = reshape(data_y_y', numel(data_y_y),1);data_y_y(isnan(data_y_y)) = [];
            end
            flag_harvest_level_x = 0;
            flag_harvest_level_y = 1;
        end
        if ~contains(bin_list(counter_bin_list), 'combined')
            if flag_harvest_level_x && flag_harvest_level_y
                bin = (count.(bin_list(counter_bin_list)));
                dir_harvest_bin = (count.dir_harvest_bin);
            elseif ~flag_harvest_level_x && ~flag_harvest_level_y
                bin = (count.([char(bin_list(counter_bin_list)) '_']));
                dir_harvest_bin = (count.(['dir_harvest_bin_']));
            elseif ~flag_harvest_level_x && flag_harvest_level_y
                bin = (count.([char(bin_list(counter_bin_list)) '_']));
                dir_harvest_bin = (count.(['dir_harvest_bin_']));
            end
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, dir_harvest_bin);
            end
            bin = cell2mat(bin);
        else
            bin = ones(length(data_x),1);
        end
        hold on
        for counter_bin = 1 : max(bin)
            if ~contains(variable_list_y(counter_variable),'var_err_eye_tgt')
                data_x_ = data_x(bin == counter_bin,:);
                data_y_ = data_y(bin == counter_bin,:);
                if ~flag_harvest_level_x && ~flag_harvest_level_y
                    data_x_ = reshape(data_x_', numel(data_x_),1);
                    data_y_ = reshape(data_y_', numel(data_y_),1);
                end
                is_nan = isnan(data_x_) | isnan(data_y_);
                data_x_(is_nan) = [];data_y_(is_nan) = [];
                bin_data_x = discretize(data_x_,edges);
                for counter_bins = 1:length(bins)
                    data_all{counter_bin, counter_bins} = data_y_(bin_data_x == counter_bins);
                    mean_data_y_bin(counter_bins) = nanmean(data_y_(bin_data_x == counter_bins)) ;
                    sem_data_y_bin(counter_bins) = nanstd(data_y_(bin_data_x == counter_bins))/sqrt(sum(bin_data_x == counter_bins)) ;
                end
            else
                data_x_ = data_x(bin == counter_bin,:);
                data_x_y_ = data_x_y(bin == counter_bin,:);
                data_y_y_ = data_y_y(bin == counter_bin,:);
                if ~flag_harvest_level_x && ~flag_harvest_level_y
                    data_x_ = reshape(data_x_', numel(data_x_),1);
                    data_x_y_ = reshape(data_x_y_', numel(data_x_y_),1);
                    data_y_y_ = reshape(data_y_y_', numel(data_y_y_),1);
                end
                is_nan = isnan(data_x_) | isnan(data_x_y_) | isnan(data_y_y_);
                data_x_(is_nan) = [];data_x_y_(is_nan) = []; data_y_y_(is_nan) = [];
                bin_data_x = discretize(data_x_,edges);
                for counter_bins = 1:length(bins)
                    [mean_data_y_bin(counter_bins), sem_data_y_bin(counter_bins),data_all{counter_bin, counter_bins}] = endpoint_error_var(data_x_y_(bin_data_x == counter_bins), data_y_y_(bin_data_x == counter_bins));
                end
            end
            %             errorbar(bins,mean_data_y_bin, sem_data_y_bin, colors(counter_bin), 'LineWidth', 2, 'MarkerSize', 10);
            shade(bins, mean_data_y_bin + sem_data_y_bin, colors(counter_bin), bins, mean_data_y_bin - sem_data_y_bin, colors(counter_bin), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(counter_bin), 'FillColor', colors(counter_bin));
            plot(bins, mean_data_y_bin, colors(counter_bin), 'LineWidth',1)

            % fit linear regression
            %             P_ = polyfit(data_x_, data_y_, 1);
            %             y_axis_hat = polyval(P_,data_x_);
            %             %plot(data_x_, data_y_,'.k')
            %             %plot(data_x_, y_axis_hat, '--r', 'LineWidth', 1)
            %             [b,~,~,~,stats] = regress(data_y_,[ones(size(data_x_)) data_x_]);
            %             r = stats(1); f = stats(2); p = stats(3); e = stats(4);
            %             test = ['r=' num2str(r), ', f=' num2str(f), ', p=' num2str(p) ', e=' num2str(e)];

            % compute correlation
            %             [r,p] = corrcoef(data_x_, data_y_);
            [r,p] = corrcoef(bins, mean_data_y_bin);
            r=r(1,2);
            p=p(1,2);
            df = length(bins)-2;
            test_r = ['r(' num2str(df) ') = ' num2str(round(r,3,'significant')), ', ' num2str(round(p,3,'significant'))];
        end

        % compute ANOVA
        %         if  ~contains(variable_list_y(counter_variable),'var_err_eye_tgt')
        [test_stat, summary, table_data] = PGH_stat_ANOVA(data_all);
        % for osf
        writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_work_Xbin_compare' '_' char(variable_list_x(counter_variable)) 'X' char(variable_list_y(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

        data_all = [];
        %         else
        %             test_stat = [];
        %             summary = [];
        %         end
        xlabel(variable_list_x(counter_variable),'Interpreter','none')
        ylabel(variable_list_y(counter_variable),'Interpreter','none')
        xlim([edges(1) edges(end)])
        ylim(y_lim)
        yticks(y_ticks)
        xticks(x_ticks)
        title_ = bin_list(counter_bin_list);
        text(max(bin)+0.1,(max(mean_data_y_bin) + min(mean_data_y_bin))/2,num2str(summary))
        if contains(bin_list(counter_bin_list), 'combined')
            test = [test_r newline test_stat];
            title([title_ test ], 'Interpreter', 'none')
        else
            test = [test_stat];
            title([title_ test ], 'Interpreter', 'none')
        end
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_work_Xbin_compare' ], 'pdf');
%close(gcf)

%%% PLOT (5) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["ITI", "eye_vigor_w", "dist_err_eye_tgt","var_err_eye_tgt", "reaction", "pupil_area"];
bin_list = ["interval_bin_", "dist_tube_bin_", "weight_bin_"];
for counter_session = 1 : count_session
    figure
    for counter_bin_list = 1 : length(bin_list)
        for counter_variable = 1 : length(variable_list)
            mean_data = [];
            sem_data = [];
            if counter_bin_list==1
                shift = 0;
            else
                shift = (counter_bin_list - 1)*length(variable_list);
            end
            y_lim = ([-inf inf]);
            y_ticks = 'auto';

            subplot(length(bin_list),length(variable_list),counter_variable + shift)
            hold on
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);
            end
            bin = cell2mat(bin);
            bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
            for counter_bin = 1 : max(bin)
                if ~contains(variable_list(counter_variable),'var_err_eye_tgt')
                    data = cell2mat(sac.(variable_list(counter_variable)));data = data(bin == counter_bin & is_bin_dir,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                    data_all{1,counter_bin} = data;
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                else
                    data_x = cell2mat(sac.err_eye_tgt_px_rot); data_x = data_x(bin == counter_bin & is_bin_dir,:); data_x = reshape(data_x', numel(data_x),1); data_x(isnan(data_x))=[];
                    data_y = cell2mat(sac.err_eye_tgt_py_rot); data_y = data_y(bin == counter_bin & is_bin_dir,:); data_y = reshape(data_y', numel(data_y),1); data_y(isnan(data_y))=[];
                    [mean_data(counter_bin), sem_data(counter_bin),data_all{1, counter_bin}] = endpoint_error_var(data_x, data_y);
                end
            end
            % compute ANOVA
            %             if  ~contains(variable_list(counter_variable),'var_err_eye_tgt')
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_sac_Xbin' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
            %             else
            %                 test = [];
            %                 summary = [];
            %             end
            errorbar(1:max(bin), mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(1), 'LineWidth',1)
            xlabel(bin_list(counter_bin_list),'Interpreter','none')
            ylabel(variable_list(counter_variable),'Interpreter', 'none')
            xlim([0.5 max(bin)+0.5])
            xticks([1:max(bin)])
            ylim(y_lim)
            yticks(y_ticks)
            title(test, 'Interpreter', 'none')
            text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
            if plot_session == 1
                % to plot session
                mean_data = [];
                sem_data = [];
                bin = count.(bin_list(counter_bin_list)){counter_session,1};
                bin_dir = count.dir_harvest_bin_{counter_session,1}; is_bin_dir = bin_dir == 1 | bin_dir == 2;
                for counter_bin = 1 : max(bin)
                    data = sac.(variable_list(counter_variable)){counter_session,1}; data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data))=[];
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                    %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
                end
                errorbar(1:max(bin), mean_data, sem_data, colors(2), 'LineWidth', 2, 'MarkerSize', 10);
                %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
                plot(1:max(bin), mean_data, colors(2), 'LineWidth',1)
                ext = [ '_' num2str(counter_session)];
                path_ext =  ['SESSION' filesep 'work' filesep ];
            end
        end
    end
    ESN_Beautify_Plot(gcf, [20,10])
    saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_sac_Xbin' ext], 'pdf');
    %close(gcf)
end
%%% PLOT (5.1) %%%
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);
        end
        [bin_dist] = disclude_data(count.dist_tube_bin_, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);

        %         % for RANOVA
        %         if ~contains(bin_list(counter_bin_list),'dist') && ~contains(variable_list(counter_variable),'var_err_eye_tgt')
        %             [data_all] = build_data_all(sac.(variable_list(counter_variable)), bin_dist ,bin);
        %             % compute ANOVA
        %             [test, summary] = PGH_stat_RANOVA(data_all);
        %             data_all = [];
        %         else
        %             test = [];
        %             summary = [];
        %         end

        bin = cell2mat(bin);
        bin_dist = cell2mat(bin_dist);
        tag_sac = cell2mat(sac.tag_sac); is_tag_nan = tag_sac == 1;
        bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
        for counter_bin_dist = 1 : max(bin_dist)
            for counter_bin = 1 : max(bin)
                if ~contains(variable_list(counter_variable),'var_err_eye_tgt')
                    data = cell2mat(sac.(variable_list(counter_variable)));data(is_tag_nan) = nan;data = data(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                    data_all{counter_bin_dist, counter_bin} = data;
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                else
                    data_x = cell2mat(sac.err_eye_tgt_px_rot); data_x = data_x(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:); data_x = reshape(data_x', numel(data_x),1); data_x(isnan(data_x))=[];
                    data_y = cell2mat(sac.err_eye_tgt_py_rot); data_y = data_y(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:); data_y = reshape(data_y', numel(data_y),1); data_y(isnan(data_y))=[];
                    if ~isempty(data_x)
                        [mean_data(counter_bin), sem_data(counter_bin),data_all{counter_bin_dist, counter_bin}] = endpoint_error_var(data_x, data_y);
                    else
                        mean_data(counter_bin) = nan;
                        sem_data(counter_bin) = nan;
                    end
                end

            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_bin_dist), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(counter_bin_dist), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_bin_dist), 'LineWidth',1)
        end
        % compute ANOVA
        if  ~contains(variable_list(counter_variable),'var_err_eye_tgt')
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_sac_Xbin_xdist' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
        else
            test = [];
            summary = [];
        end
        xlabel(bin_list(counter_bin_list),'Interpreter','none');
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
        title(test, 'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_sac_Xbin_xdist' ext], 'pdf');
%close(gcf)

%%% PLOT (6) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["ITI", "eye_vigor_w", "dist_err_eye_tgt","var_err_eye_tgt", "reaction", "pupil_area"];
bin_list = ["combined", "interval_bin_", "dist_tube_bin_", "weight_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        title_ = bin_list(counter_bin_list);
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        if ~contains(bin_list(counter_bin_list), 'combined')
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);
            end
            % for RANOVA
            if ~contains(variable_list(counter_variable),'var_err_eye_tgt') && ~contains(bin_list(counter_bin_list), 'combined')
                [data_all] = build_data_all(sac.(variable_list(counter_variable)), [] ,bin);
                % compute ANOVA
                [test, summary, table_data] = PGH_stat_RANOVA(data_all);
                % for osf
                writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_sac_Xsac' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

                data_all = [];
            else
                test = [];
                summary = [];
            end
            bin = cell2mat(bin);

        else
            bin = ones(length(cell2mat(count.(bin_list(2)))),1);
        end

        bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
        for counter_bin = 1 : max(bin)
            if ~contains(variable_list(counter_variable),'var_err_eye_tgt')
                data = cell2mat(sac.(variable_list(counter_variable)));data = data(bin == counter_bin & is_bin_dir,:);
                % for anova
                if contains(bin_list(counter_bin_list), 'combined')
                    for counter_event = 1 : size(data,2)
                        data_all{counter_bin, counter_event} = data(:,counter_event);
                    end
                end
                mean_data = (nanmean(data)) ;
                sem_data = nanstd(data)/sqrt(size(data,1)) ;
            else
                data_x = cell2mat(sac.err_eye_tgt_px_rot); data_x = data_x(bin == counter_bin & is_bin_dir,:);
                data_y = cell2mat(sac.err_eye_tgt_py_rot); data_y = data_y(bin == counter_bin & is_bin_dir,:);
                if ~isempty(data_x)
                    [mean_data, sem_data] = endpoint_error_var(data_x, data_y);
                else
                    mean_data = nan(1,size(data_x,2));
                    sem_data = nan(1,size(data_x,2));
                end
            end
            %             errorbar(mean_data, sem_data, colors(counter_bin), 'LineWidth', 2, 'MarkerSize', 10);
            shade(1:length(mean_data), mean_data + sem_data, colors(counter_bin), 1:length(mean_data), mean_data - sem_data, colors(counter_bin), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(counter_bin), 'FillColor', colors(counter_bin));
            plot(1:length(mean_data), mean_data, colors(counter_bin), 'LineWidth',1)
        end
        % compute ANOVA
        if  ~contains(variable_list(counter_variable),'var_err_eye_tgt') && contains(bin_list(counter_bin_list),'combined')
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_sac_Xsac' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
        end
        %         summary = [];
        xlim([0.5 10.5])
        xticks([1:10])
        xlabel('num sac in work')
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        ylim(y_lim)
        yticks(y_ticks)
        title_ = [title_ test];
        title(title_,'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_sac_Xsac' ], 'pdf');
%close(gcf)

%%% PLOT (7) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["ITI", "eye_vigor_w", "dist_err_eye_tgt","var_err_eye_tgt", "reaction", "pupil_area"];
bin_list = ["interval_bin_", "dist_tube_bin_", "weight_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);
        end
        bin = cell2mat(bin);
        for counter_tag_sac = 1 : count.num_tag_sac-6
            for counter_bin = 1 : max(bin)
                tag_sac = cell2mat(sac.tag_sac);
                if counter_tag_sac == 1
                    is_tag_sac = tag_sac == 1;
                elseif counter_tag_sac == 2
                    is_tag_sac = tag_sac == 4 ;
                elseif counter_tag_sac == 3
                    is_tag_sac = tag_sac == 6 ;
                elseif counter_tag_sac == 4
                    is_tag_sac = tag_sac == 10;
                end
                if ~contains(variable_list(counter_variable),'var_err_eye_tgt')
                    data = cell2mat(sac.(variable_list(counter_variable)));data(~is_tag_sac) = nan; data = data(bin == counter_bin,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                else
                    data_x = cell2mat(sac.err_eye_tgt_px_rot); data_x(~is_tag_sac) = nan; data_x = data_x(bin == counter_bin,:); data_x = reshape(data_x', numel(data_x),1); data_x(isnan(data_x))=[];
                    data_y = cell2mat(sac.err_eye_tgt_py_rot); data_y(~is_tag_sac) = nan; data_y = data_y(bin == counter_bin,:); data_y = reshape(data_y', numel(data_y),1); data_y(isnan(data_y))=[];
                    if ~isempty(data_x)
                        [mean_data(counter_bin), sem_data(counter_bin)] = endpoint_error_var(data_x, data_y);
                    else
                        mean_data(counter_bin) = nan;
                        sem_data(counter_bin) = nan;
                    end
                end
            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_tag_sac), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_tag_sac), 'LineWidth',1)
        end
        xlabel(bin_list(counter_bin_list),'Interpreter','none')
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_sac_Xbin_type' ], 'pdf');
%close(gcf)

%%% PLOT (8) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['WORK' filesep ];
variable_list = ["ITI", "eye_vigor_w", "dist_err_eye_tgt","var_err_eye_tgt", "reaction", "pupil_area"];
bin_list = ["combined", "interval_bin_", "dist_tube_bin_", "weight_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        title_ = bin_list(counter_bin_list);
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        if ~contains(bin_list(counter_bin_list), 'combined')
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_sac_l, count.disclude_sac_r, count.dir_harvest_bin_);
            end
            bin = cell2mat(bin);
        else
            bin = ones(length(cell2mat(count.(bin_list(2)))),1);
        end
        for counter_tag_sac = 1 : count.num_tag_sac-6
            for counter_bin = 1 : max(bin)
                tag_sac = cell2mat(sac.tag_sac);
                if counter_tag_sac == 1
                    is_tag_sac = tag_sac == 1;
                elseif counter_tag_sac == 2
                    is_tag_sac = tag_sac == 4 ;
                elseif counter_tag_sac == 3
                    is_tag_sac = tag_sac == 6 ;
                elseif counter_tag_sac == 4
                    is_tag_sac = tag_sac == 10;
                end
                if ~contains(variable_list(counter_variable),'var_err_eye_tgt')
                    data = cell2mat(sac.(variable_list(counter_variable)));data(~is_tag_sac) = nan; data = data(bin == counter_bin,:);
                    mean_data = (nanmean(data)) ;
                    sem_data = nanstd(data)/sqrt(length(data)) ;
                else
                    data_x = cell2mat(sac.err_eye_tgt_px_rot); data_x(~is_tag_sac) = nan; data_x = data_x(bin == counter_bin,:);
                    data_y = cell2mat(sac.err_eye_tgt_py_rot); data_y(~is_tag_sac) = nan; data_y = data_y(bin == counter_bin,:);
                    [mean_data, sem_data] = endpoint_error_var(data_x, data_y);
                end

            end
            errorbar(1:length(mean_data), mean_data, sem_data, colors(counter_tag_sac), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:length(mean_data), mean_data, colors(counter_tag_sac), 'LineWidth',1)
        end
        xlim([0.5 10.5])
        xticks([1:10])
        xlabel('num sac in work')
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        ylim(y_lim)
        yticks(y_ticks)
        title(title_,'Interpreter', 'none')
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_sac_Xsac_type' ], 'pdf');
%close(gcf)


%% Figure - harvest period
% %%% PLOT (1)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
figure
hold on
for counter_session = 1 : count.num_session
    mode_data(counter_session,1) = mode(count.dist_tube_bin{counter_session, 1});
end
x_axis = categorical({'near', 'mid', 'far'}); x_axis = reordercats(x_axis,{'near', 'mid', 'far'});
data = count.num_lick;
for counter_bin = 1 : max(mode_data)
    data_{counter_bin} = data(mode_data==counter_bin);
    mean_data(counter_bin) = nanmean(data_{counter_bin});
    sem_data(counter_bin) = nanstd(data_{counter_bin})/sqrt(length(data_{counter_bin}));
    %     plot(x_axis(counter_bin),data_{counter_bin},'or','MarkerSize', 5)
end
errorbar(x_axis, mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
plot(x_axis, mean_data, colors(1), 'LineWidth',1)
xlabel('tube dist')
ylabel('# licks in session')
title([num2str(mean_data(1)) '+/-' num2str(sem_data(1)) ' | ' num2str(mean_data(2)) '+/-' num2str(sem_data(2)) ' | ' num2str(mean_data(3)) '+/-' num2str(sem_data(3))])
ESN_Beautify_Plot(gcf, [5,5])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_licks' ], 'pdf');
%close(gcf)

%%% PLOT (2)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["num_lick_harvest", "duration_harvest", "rew_str_harvest", "rew_end_harvest", "rew_consumed_harvest"];
figure
for counter_variable = 1 : length(variable_list)
    data = cell2mat(lick.(variable_list(counter_variable)));data = reshape(data', numel(data),1);
     % for osf           
    writetable(array2table(data), [path.out_path 'data' filesep path.path_data_monkey_sorted '_harvest_hist_' char(variable_list(counter_variable))  ext '.csv']);

    mean_data = nanmean(data) ;
    std_data = nanstd(data)/sqrt(sum(~isnan(data)));
    if contains(variable_list(counter_variable), 'num_lick_harvest')
        x_lim = [0 80];
        x_ticks = [0:10:80];
        data(data>80) = nan;
    elseif contains(variable_list(counter_variable), 'duration_harvest')
        x_lim = [0 60];
        x_ticks = [0:10:60];
        data(data>60) = nan;
    elseif contains(variable_list(counter_variable), 'rew_str_harvest')
        x_lim = [0 1.25];
        x_ticks = [0:0.25:1.25];
    elseif contains(variable_list(counter_variable), 'rew_end_harvest')
        x_lim = [0 1.25];
        x_ticks = [0:0.25:1.25];
    elseif contains(variable_list(counter_variable), 'rew_consumed_harvest')
        x_lim = [0 1.25];
        x_ticks = [0:0.25:1.25];
    end
    subplot(1,length(variable_list),counter_variable)
    hold on
    histogram(data,'FaceColor', 'k','Normalization','probability')
    xline(mean_data,'-r', 'LineWidth',2)
    xline(mean_data+std_data,'--r', 'LineWidth',1)
    xline(mean_data-std_data,'--r', 'LineWidth',1)
    xlabel(variable_list(counter_variable),'Interpreter', 'none')
    ylabel('prob.')
    xlim(x_lim)
    xticks(x_ticks)
    title([num2str(mean_data) '+/-' num2str(std_data)])
end
sgtitle(['Total trials: ' num2str(nanmean(count.num_trial)) ' +/- ' num2str(nanstd(count.num_trial)/sqrt(length(count.num_trial))) ' | Total intervals: ' num2str(nanmean(count.num_harvest))  ' +/- ' num2str(nanstd(count.num_harvest)/sqrt(length(count.num_harvest)))])
ESN_Beautify_Plot(gcf, [20,3])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_harvest_hist' ], 'pdf');
%close(gcf)

%%% PLOT (3)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["num_lick_harvest", "duration_harvest", "rew_str_harvest", "rew_end_harvest", "rew_consumed_harvest"];
bin_list = ["interval_bin", "dist_tube_bin", "num_trial_bin", "weight_bin"];
for counter_session = 1 : count_session
    figure
    for counter_bin_list = 1 : length(bin_list)
        for counter_variable = 1 : length(variable_list)
            mean_data = [];
            sem_data = [];
            % handle shift in subplots
            if counter_bin_list==1
                shift = 0;
            else
                shift = (counter_bin_list - 1)*length(variable_list);
            end
            y_lim = ([-inf inf]);
            y_ticks = 'auto';
            subplot(length(bin_list),length(variable_list),counter_variable + shift)
            hold on
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_harvest_l, count.disclude_harvest_r, count.dir_harvest_bin);
            end
            bin = cell2mat(bin);
            bin_dir = cell2mat(count.dir_harvest_bin); is_bin_dir = bin_dir == 1 | bin_dir == 2;
            for counter_bin = 1 : max(bin)
                data = cell2mat(lick.(variable_list(counter_variable)));data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data))=[];
                data_all{1,counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
            end
            % compute ANOVA
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);

            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_harvest_Xbin' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
            errorbar(1:max(bin), mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(1), 'LineWidth',1)
            xlabel(bin_list(counter_bin_list),'Interpreter','none');
            ylabel(variable_list(counter_variable),'Interpreter', 'none')
            xlim([0.5 max(bin)+0.5])
            xticks([1:max(bin)])
            ylim(y_lim)
            yticks(y_ticks)
            title(test, 'Interpreter', 'none')
            text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
            if plot_session == 1
                % to plot session
                mean_data = [];
                sem_data = [];
                bin = count.(bin_list(counter_bin_list)){counter_session,1};
                bin_dir = count.dir_harvest_bin{counter_session,1}; is_bin_dir = bin_dir == 1 | bin_dir == 2;
                for counter_bin = 1 : max(bin)
                    data = lick.(variable_list(counter_variable)){counter_session,1}; data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data) | data<=0)=[];
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                    %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
                end
                errorbar(1:max(bin), mean_data, sem_data, colors(2), 'LineWidth', 2, 'MarkerSize', 10);
                %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
                plot(1:max(bin), mean_data, colors(2), 'LineWidth',1)
                ext = [ '_' num2str(counter_session)];
                path_ext =  ['SESSION' filesep 'harvest' filesep ];
            end
        end
    end
    ESN_Beautify_Plot(gcf, [20,10])
    saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_harvest_Xbin' ext], 'pdf');
    %close(gcf)
end
%%% PLOT (3.1) %%%
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_harvest_l, count.disclude_harvest_r, count.dir_harvest_bin_);
        end
        [bin_dist] = disclude_data(count.dist_tube_bin, count.disclude_harvest_l, count.disclude_harvest_r, count.dir_harvest_bin);

        %          % for RANOVA
        %         [data_all] = build_data_all(lick.(variable_list(counter_variable)), bin_dist ,bin);
        %         % compute ANOVA
        %         if ~contains(bin_list(counter_bin_list),'dist')
        %             [test, summary] = PGH_stat_RANOVA(data_all);
        %             data_all = [];
        %         else
        %             test = [];
        %             summary = [];
        %         end

        bin = cell2mat(bin);
        bin_dist = cell2mat(bin_dist);
        bin_dir = cell2mat(count.dir_harvest_bin); is_bin_dir = bin_dir == 1 | bin_dir == 2;
        for counter_bin_dist = 1 : max(bin_dist)
            for counter_bin = 1 : max(bin)
                data = cell2mat(lick.(variable_list(counter_variable))); data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:); data(isnan(data))=[];
                data_all{counter_bin_dist, counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_bin_dist), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(counter_bin_dist), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_bin_dist), 'LineWidth',1)
        end
        % compute ANOVA
        [test, summary, table_data] = PGH_stat_ANOVA(data_all);
        % for osf
        writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_harvest_Xbin_xdist' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

        data_all = [];
        xlabel(bin_list(counter_bin_list),'Interpreter','none');
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
        title(test, 'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_harvest_Xbin_xdist' ext], 'pdf');
%close(gcf)

%%% PLOT (4)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list_x = [ "pupil_area",  "pupil_area"];
variable_list_y = [ "tongue_vigor_pro", "duration_harvest"];
bin_list = ["combined", "interval_bin","dist_tube_bin", "num_trial_bin", "weight_bin"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list_x)
        mean_data_y_bin = [];
        sem_data_y_bin = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list_x);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        if contains(variable_list_x(counter_variable), 'num_lick_harvest')
            bins = 3 : 30;
            edges = 2.5 : 1 : 30.5;
            x_ticks = [1 10 20 30];
        elseif contains(variable_list_x(counter_variable), 'rew_str_harvest')
            bins = 0 :0.1: 1.25;
            edges = 0 : 0.1 : 1.25;
            x_ticks = [0 0.5 1 1.25];
        elseif contains(variable_list_x(counter_variable), 'pupil_area') %&& contains(bin_list(counter_bin_list), 'combined')
            %             if contains(path.path_data_monkey_sorted,'59d')
            %                 bins = 5750 :250: 8250;
            %                 edges = 5500 : 250 : 8500;
            %                 x_ticks = [bins];
            %             elseif contains(path.path_data_monkey_sorted,'125d')
            %                 bins = 5750 :250: 7500;
            %                 edges = 5500 : 250 : 7750;
            %                 x_ticks = [bins];
            %             end
            bins = -1 :0.2: 1;
            edges = -1.1 : 0.2 : 1.1;
            x_ticks = [bins];
        elseif contains(variable_list_x(counter_variable), 'tongue_vigor_pro') %&& contains(bin_list(counter_bin_list), 'combined')
            bins = 0.6 :0.1: 1.4;
            edges = 0.5 : 0.1 : 1.5;
            x_ticks = [bins];
        end
        subplot(length(bin_list) ,length(variable_list_x),counter_variable + shift)
        data_x = cell2mat(lick.(variable_list_x(counter_variable)));
        data_y = cell2mat(lick.(variable_list_y(counter_variable)));
        if size(data_x,1) == count.num_session && size(data_y,1) == count.num_session % both harvest level
            data_x = reshape(data_x', numel(data_x),1);
            data_y = reshape(data_y', numel(data_y),1);
            flag_harvest_level_x = 1;
            flag_harvest_level_y = 1;
        elseif size(data_x,1) ~= count.num_session && size(data_y,1) ~= count.num_session % both lick level
            flag_harvest_level_x = 0;
            flag_harvest_level_y = 0;
        elseif size(data_x,1) ~= count.num_session && size(data_y,1) == count.num_session % x lick level, y harvest  level
            data_x = nanmean(data_x,2);
            data_y = reshape(data_y', numel(data_y),1);data_y(isnan(data_y)) = [];
            flag_harvest_level_x = 0;
            flag_harvest_level_y = 1;
        end
        if ~contains(bin_list(counter_bin_list), 'combined')
            if flag_harvest_level_x && flag_harvest_level_y
                bin = (count.(bin_list(counter_bin_list)));
                dir_harvest_bin = (count.dir_harvest_bin);
            elseif ~flag_harvest_level_x && ~flag_harvest_level_y
                bin = (count.([char(bin_list(counter_bin_list)) '_']));
                dir_harvest_bin = (count.(['dir_harvest_bin_']));
            elseif ~flag_harvest_level_x && flag_harvest_level_y
                bin = (count.([char(bin_list(counter_bin_list)) '_']));
                dir_harvest_bin = (count.(['dir_harvest_bin_']));
            end
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, dir_harvest_bin);
            end
            bin = cell2mat(bin);
        else
            bin = ones(length(data_x),1);
        end

        hold on
        for counter_bin = 1 : max(bin)
            data_x_ = data_x(bin == counter_bin,:);
            data_y_ = data_y(bin == counter_bin,:);
            if ~flag_harvest_level_x && ~flag_harvest_level_y
                data_x_ = reshape(data_x_', numel(data_x_),1);
                data_y_ = reshape(data_y_', numel(data_y_),1);
            end
            is_nan = isnan(data_x_) | isnan(data_y_);
            data_x_(is_nan) = [];data_y_(is_nan) = [];
            bin_data_x = discretize(data_x_,edges);
            for counter_bins = 1:length(bins)
                data_all{counter_bin, counter_bins} = data_y_(bin_data_x == counter_bins);
                mean_data_y_bin(counter_bins) = nanmean(data_y_(bin_data_x == counter_bins)) ;
                sem_data_y_bin(counter_bins) = nanstd(data_y_(bin_data_x == counter_bins))/sqrt(sum(bin_data_x == counter_bins)) ;
            end

            %             errorbar(bins,mean_data_y_bin, sem_data_y_bin, colors(counter_bin), 'LineWidth', 2, 'MarkerSize', 10);
            shade(bins, mean_data_y_bin + sem_data_y_bin, colors(counter_bin), bins, mean_data_y_bin - sem_data_y_bin, colors(counter_bin), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(counter_bin), 'FillColor', colors(counter_bin));
            plot(bins, mean_data_y_bin, colors(counter_bin), 'LineWidth',1)
            % fit linear regression
            %             P_ = polyfit(data_x_, data_y_, 1);
            %             y_axis_hat = polyval(P_,data_x_);
            %             %plot(data_x_, data_y_,'.k')
            %             %plot(data_x_, y_axis_hat, '--r', 'LineWidth', 1)
            %             [b,~,~,~,stats] = regress(data_y_,[ones(size(data_x_)) data_x_]);
            %             r = stats(1); f = stats(2); p = stats(3); e = stats(4);
            %             test = ['r=' num2str(r), ', f=' num2str(f), ', p=' num2str(p) ', e=' num2str(e)];

            % compute correlation
            %[r,p] = corrcoef(data_x_, data_y_);
            [r,p] = corrcoef(bins, mean_data_y_bin);
            r=r(1,2);
            p=p(1,2);
            df = length(bins)-2;
            test_r = ['r(' num2str(df) ') = ' num2str(r), ', p=' num2str(p)];
        end
        % compute RANOVA
        [test_stat, summary, table_data] = PGH_stat_ANOVA(data_all);
        % for osf
        writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_harvest_Xbin_xdist' '_' char(variable_list_x(counter_variable)) 'X' char(variable_list_y(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

        data_all = [];

        xlabel(variable_list_x(counter_variable),'Interpreter','none')
        ylabel(variable_list_y(counter_variable),'Interpreter','none')
        xlim([edges(1) edges(end)])
        xticks(x_ticks)
        ylim(y_lim)
        yticks(y_ticks)
        title_ = bin_list(counter_bin_list);
        text(max(bin)+0.1,(max(mean_data_y_bin) + min(mean_data_y_bin))/2,num2str(summary))
        if contains(bin_list(counter_bin_list), 'combined')
            test = [test_r newline test_stat];
            title([title_ test ], 'Interpreter', 'none')
        else
            test = [test_stat];
            title([title_ test ], 'Interpreter', 'none')
        end
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_harvest_Xbin_compare' ], 'pdf');
%close(gcf)

%%% PLOT (5) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["ILI", "tongue_vigor_pro", "tongue_vigor_ret","dist_err_tongue_rew", "tongue_duration", "rew", "pupil_area", "J"];
bin_list = ["interval_bin_", "dist_tube_bin_", "num_trial_bin_", "weight_bin_"];
for counter_session = 1 : count_session
    figure
    for counter_bin_list = 1 : length(bin_list)
        for counter_variable = 1 : length(variable_list)
            mean_data = [];
            sem_data = [];
            if counter_bin_list==1
                shift = 0;
            else
                shift = (counter_bin_list - 1)*length(variable_list);
            end
            y_lim = ([-inf inf]);
            y_ticks = 'auto';
            subplot(length(bin_list),length(variable_list),counter_variable + shift)
            hold on
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
            end
            bin = cell2mat(bin);
            tag_lick = cell2mat(lick.tag_lick); is_tag_nan = tag_lick == 1;
            bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
            for counter_bin = 1 : max(bin)
                data = cell2mat(lick.(variable_list(counter_variable)));
                data(is_tag_nan) = nan;
                data = data(bin == counter_bin & is_bin_dir,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                data_all{1,counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
            end
            % compute ANOVA
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_lick_Xbin' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
            %             if (contains(bin_list(counter_bin_list),'weight') || contains(bin_list(counter_bin_list),'num_trial') || contains(bin_list(counter_bin_list),'dist'))  && (contains(variable_list(counter_variable),'pupil_area') || contains(variable_list(counter_variable),'tongue_vigor_pro'))
            %                 data = cell2mat(lick.(variable_list(counter_variable)));
            %                 data_1 = data(bin == 1 & is_bin_dir,:);data_1 = reshape(data_1', numel(data_1),1); data_1(isnan(data_1))=[];
            %                 data_2 = data(bin == 2 & is_bin_dir,:);data_2 = reshape(data_2', numel(data_2),1); data_2(isnan(data_2))=[];
            %                 data_3 = data(bin == 3 & is_bin_dir,:);data_3 = reshape(data_3', numel(data_3),1); data_3(isnan(data_3))=[];
            %                 data_4 = data(bin == 4 & is_bin_dir,:);data_4 = reshape(data_4', numel(data_4),1); data_4(isnan(data_4))=[];
            %                 data_5 = data(bin == 5 & is_bin_dir,:);data_5 = reshape(data_5', numel(data_5),1); data_5(isnan(data_5))=[];
            %
            %                 mean_data_1 = nanmean(data_1);
            %                 CI_data_1 = nanstd(data_1)/sqrt(length(data_1))*1.96;
            %                 mean_data_2 = nanmean(data_2);
            %                 CI_data_2 = nanstd(data_2)/sqrt(length(data_2))*1.96;
            %                 mean_data_3 = nanmean(data_3);
            %                 CI_data_3 = nanstd(data_3)/sqrt(length(data_3))*1.96;
            %                 mean_data_4 = nanmean(data_4);
            %                 CI_data_4 = nanstd(data_4)/sqrt(length(data_4))*1.96;
            %                 mean_data_5 = nanmean(data_5);
            %                 CI_data_5 = nanstd(data_5)/sqrt(length(data_5))*1.96;
            %
            %                 test = [num2str(mean_data_1) '+/-' num2str(CI_data_1) newline ...
            %                     num2str(mean_data_2) '+/-' num2str(CI_data_2) newline ...
            %                     num2str(mean_data_3) '+/-' num2str(CI_data_3) newline ...
            %                     num2str(mean_data_4) '+/-' num2str(CI_data_4) newline ...
            %                     num2str(mean_data_5) '+/-' num2str(CI_data_5) newline];
            %             else
            %                 test = [];
            %             end
            errorbar(1:max(bin), mean_data, sem_data, colors(1), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(1), 'LineWidth',1)
            xlabel(bin_list(counter_bin_list),'Interpreter','none');
            ylabel(variable_list(counter_variable),'Interpreter', 'none')
            xlim([0.5 max(bin)+0.5])
            xticks([1:max(bin)])
            ylim(y_lim)
            yticks(y_ticks)
            title(test, 'Interpreter', 'none')
            text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
            if plot_session == 1
                % to plot session
                mean_data = [];
                sem_data = [];
                bin = count.(bin_list(counter_bin_list)){counter_session,1};
                bin_dir = count.dir_harvest_bin_{counter_session,1}; is_bin_dir = bin_dir == 1 | bin_dir == 2;
                for counter_bin = 1 : max(bin)
                    data = lick.(variable_list(counter_variable)){counter_session,1}; data = reshape(data', numel(data),1); data = data(bin == counter_bin & is_bin_dir,:); data(isnan(data))=[];
                    mean_data(counter_bin) = (nanmean(data)) ;
                    sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
                    %plot(counter_bin, data, ['o' colors(2)], 'MarkerSize', 5)
                end
                errorbar(1:max(bin), mean_data, sem_data, colors(2), 'LineWidth', 2, 'MarkerSize', 10);
                %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
                plot(1:max(bin), mean_data, colors(2), 'LineWidth',1)
                ext = [ '_' num2str(counter_session)];
                path_ext =  ['SESSION' filesep 'harvest' filesep ];
            end
        end
    end
    ESN_Beautify_Plot(gcf, [20,10])
    saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_Xbin' ext], 'pdf');
    %close(gcf)
end
%%% PLOT (5.1) %%%
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
        end
        [bin_dist] = disclude_data(count.dist_tube_bin_, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);

        % for RANOVA
        %         if ~contains(bin_list(counter_bin_list),'dist')
        %             [data_all] = build_data_all(lick.(variable_list(counter_variable)), bin_dist ,bin);
        %             % compute ANOVA
        %             [test, summary] = PGH_stat_RANOVA(data_all);
        %             data_all = [];
        %         else
        %             test = [];
        %             summary = [];
        %         end
        bin = cell2mat(bin);
        bin_dist = cell2mat(bin_dist);
        tag_lick = cell2mat(lick.tag_lick); is_tag_nan = tag_lick == 1;
        bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
        for counter_bin_dist = 1 : max(bin_dist)
            for counter_bin = 1 : max(bin)
                data = cell2mat(lick.(variable_list(counter_variable)));
                data(is_tag_nan) = nan;
                data = data(bin == counter_bin & is_bin_dir & bin_dist == counter_bin_dist,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                data_all{counter_bin_dist, counter_bin} = data;
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_bin_dist), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(counter_bin_dist), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_bin_dist), 'LineWidth',1)
        end
        % compute ANOVA
        [test, summary, table_data] = PGH_stat_ANOVA(data_all);
        % for osf
        writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_lick_Xbin_xdist' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

        data_all = [];
        xlabel(bin_list(counter_bin_list),'Interpreter','none');
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
        title(test, 'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_Xbin_xdist' ext], 'pdf');
%close(gcf)

%%% PLOT (6) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["ILI", "tongue_vigor_pro", "tongue_vigor_ret","dist_err_tongue_rew", "tongue_duration", "rew", "pupil_area"];
bin_list = ["combined", "interval_bin_", "dist_tube_bin_", "num_trial_bin_", "weight_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        title_ = bin_list(counter_bin_list);
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        if ~contains(bin_list(counter_bin_list), 'combined')
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
            end
            % for RANOVA
            if ~contains(bin_list(counter_bin_list), 'combined')
                [data_all] = build_data_all(lick.(variable_list(counter_variable)), [] ,bin);
                % compute ANOVA
                [test, summary, table_data] = PGH_stat_RANOVA(data_all);
                % for osf
                writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_lick_Xlick' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

                data_all = [];
            else
                test = [];
                summary = [];
            end
            bin = cell2mat(bin);
        else
            bin = ones(length(cell2mat(count.(bin_list(2)))),1);
        end
        tag_lick = cell2mat(lick.tag_lick); is_tag_nan = tag_lick == 1;
        bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
        for counter_bin = 1 : max(bin)
            data = cell2mat(lick.(variable_list(counter_variable)));
            data(is_tag_nan) = nan;
            data = data(bin == counter_bin & is_bin_dir,:);
            for counter_event = 1 : size(data,2)
                data_all{counter_bin, counter_event} = data(:,counter_event);
            end
            mean_data = (nanmean(data)) ;
            sem_data = nanstd(data)/sqrt(size(data,1)) ;
            %             errorbar(mean_data, sem_data, colors(counter_bin), 'LineWidth', 2, 'MarkerSize', 10);
            shade(1:length(mean_data), mean_data + sem_data, colors(counter_bin), 1:length(mean_data), mean_data - sem_data, colors(counter_bin), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(counter_bin), 'FillColor', colors(counter_bin));
            plot(1:length(mean_data), mean_data, colors(counter_bin), 'LineWidth',1)
        end
        % compute RANOVA
        if contains(bin_list(counter_bin_list),'combined')
            [test, summary, table_data] = PGH_stat_ANOVA(data_all);
            % for osf
            writetable(table_data, [path.out_path 'data' filesep path.path_data_monkey_sorted '_lick_Xlick' '_' char(variable_list(counter_variable)) '_' char(bin_list(counter_bin_list)) ext '.csv']);

            data_all = [];
        end
        summary = [];
        xlim([0.5 30.5])
        xticks([1 10 20 30])
        xlabel('num lick in harvest')
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        ylim(y_lim)
        title_ = [title_ test];
        title(title_,'Interpreter', 'none')
        text(max(bin)+0.1,(max(mean_data) + min(mean_data))/2,num2str(summary))
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_Xlick' ], 'pdf');
%close(gcf)

%%% PLOT (6.1) MODEL %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["J"];
bin_dist_list = ["dist_tube_bin_"];
bin_trial_list = ["num_trial_bin_"];
bin_weight_list = ["weight_bin_"];
figure
if contains(bin_dist_list,'dist') && plot_session == 0
    [bin_dist] = disclude_data(count.(bin_dist_list), count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
end
bin_dist = cell2mat(bin_dist);
bin_trial = cell2mat(count.(bin_trial_list));
bin_weight =  cell2mat(count.(bin_weight_list));
tag_lick = cell2mat(lick.tag_lick); is_tag_nan = tag_lick == 1;
bin_dir = cell2mat(count.dir_harvest_bin_); is_bin_dir = bin_dir == 1 | bin_dir == 2;
for counter_variable = 1 : length(variable_list)
    for counter_weight = 1 : max(bin_weight)
        mean_data = [];
        sem_data = [];
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        if counter_weight == 1
            shift = 0;
        else
            shift = max(bin_dist);
        end
        for counter_dist = 1 : max(bin_dist)
            for counter_trial = 1 : max(bin_trial)
                title_ = [char(bin_dist_list) ': ' num2str(counter_dist) newline char(bin_weight_list) ': ' num2str(counter_weight)];
                subplot(max(bin_weight),max(bin_dist),counter_dist + shift)
                hold on
                data = cell2mat(lick.(variable_list(counter_variable)));
                data(is_tag_nan) = nan;
                data = data(bin_dist == counter_dist & bin_weight == counter_weight & bin_trial == counter_trial & is_bin_dir,:);
                mean_data = (nanmean(data)) ;
                sem_data = nanstd(data)/sqrt(size(data,1)) ;
                %             errorbar(mean_data, sem_data, colors(counter_bin), 'LineWidth', 2, 'MarkerSize', 10);
                shade(1:length(mean_data), mean_data + sem_data, colors(counter_trial), 1:length(mean_data), mean_data - sem_data, colors(counter_trial), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(counter_trial), 'FillColor', colors(counter_trial));
                plot(1:length(mean_data), mean_data, colors(counter_trial), 'LineWidth',1)
            end
            xlim([0.5 30.5])
            xticks([1 10 20 30])
            xlabel('num lick in harvest')
            ylabel(variable_list(counter_variable),'Interpreter', 'none')
            ylim(y_lim)
            title(title_,'Interpreter', 'none')
        end
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_Xlick_MODEL' ], 'pdf');
%close(gcf)

%%% PLOT (7) %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = ["ILI", "tongue_vigor_pro", "tongue_vigor_ret","dist_err_tongue_rew", "tongue_duration", "rew", "pupil_area"];
bin_list = ["interval_bin_", "dist_tube_bin_", "num_trial_bin_", "weight_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';
        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        bin = (count.(bin_list(counter_bin_list)));
        if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
            [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
        end
        bin = cell2mat(bin);
        for counter_tag_lick = 1 : count.num_tag_lick
            for counter_bin = 1 : max(bin)
                data = cell2mat(lick.(variable_list(counter_variable)));
                tag_lick = cell2mat(lick.tag_lick);
                if counter_tag_lick == 1
                    is_tag_lick = tag_lick == 1;
                elseif counter_tag_lick == 2
                    is_tag_lick = tag_lick == 2 | tag_lick == 3;
                elseif counter_tag_lick == 4
                    is_tag_lick = tag_lick == 4 | tag_lick == 5;
                elseif counter_tag_lick == 6
                    is_tag_lick = tag_lick == 6 | tag_lick == 7;
                else
                    is_tag_lick = tag_lick == 0;
                end
                data(~is_tag_lick) = nan;
                data = data(bin == counter_bin,:);data = reshape(data', numel(data),1); data(isnan(data))=[];
                mean_data(counter_bin) = (nanmean(data)) ;
                sem_data(counter_bin) = nanstd(data)/sqrt(length(data)) ;
            end
            errorbar(1:max(bin), mean_data, sem_data, colors(counter_tag_lick), 'LineWidth', 2, 'MarkerSize', 10);
            %shade(1:max(bin), mean_data + sem_data, colors(1), 1:max(bin), mean_data - sem_data, colors(1), 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',colors(1), 'FillColor', colors(1));
            plot(1:max(bin), mean_data, colors(counter_tag_lick), 'LineWidth',1)
        end
        xlabel(bin_list(counter_bin_list),'Interpreter','none');
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        xlim([0.5 max(bin)+0.5])
        xticks([1:max(bin)])
        ylim(y_lim)
        yticks(y_ticks)
    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_Xbin_type' ], 'pdf');
%close(gcf)
%% Figure - RPE
%%% PLOT (1)  %%%
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];
variable_list = [ "tongue_vigor_pro", "pupil_area"];
bin_list = ["combined", "interval_bin_", "dist_tube_bin_", "num_trial_bin_"];
figure
for counter_bin_list = 1 : length(bin_list)
    for counter_variable = 1 : length(variable_list)
        mean_data = [];
        sem_data = [];
        if counter_bin_list==1
            shift = 0;
        else
            shift = (counter_bin_list - 1)*length(variable_list);
        end
        y_lim = ([-inf inf]);
        y_ticks = 'auto';

        x_label = categorical({'T prev', 'T curr', 'F prev' 'F curr', 'T diff', 'F diff'});
        x_label = reordercats(x_label,{'T prev', 'T curr', 'F prev' 'F curr','T diff', 'F diff'});

        subplot(length(bin_list),length(variable_list),counter_variable + shift)
        hold on
        if ~contains(bin_list(counter_bin_list), 'combined')
            bin = (count.(bin_list(counter_bin_list)));
            if contains(bin_list(counter_bin_list),'dist') && plot_session == 0
                [bin] = disclude_data(bin, count.disclude_lick_l, count.disclude_lick_r, count.dir_harvest_bin_);
            end
            bin = cell2mat(bin);
        else
            bin = ones(length(cell2mat(count.(bin_list(2)))),1);
        end
        for counter_bin = 1 : max(bin)
            data = cell2mat(lick.(variable_list(counter_variable))); data = data(bin == counter_bin,:);data = reshape(data', numel(data),1);
            tag_lick = cell2mat(lick.tag_lick); tag_lick = tag_lick(bin == counter_bin,:); tag_lick = reshape(tag_lick', numel(tag_lick),1);

            %RPE
            is_TT_prev = cell2mat(count.is_TT); is_TT_prev = is_TT_prev(bin == counter_bin,1:30);is_TT_prev = reshape(is_TT_prev', numel(is_TT_prev),1);
            is_TT_curr = cell2mat(count.is_TT); is_TT_curr = is_TT_curr(bin == counter_bin,31:60);is_TT_curr = reshape(is_TT_curr', numel(is_TT_curr),1);
            is_FF_prev = cell2mat(count.is_FF); is_FF_prev = is_FF_prev(bin == counter_bin,1:30);is_FF_prev = reshape(is_FF_prev', numel(is_FF_prev),1);
            is_FF_curr = cell2mat(count.is_FF); is_FF_curr = is_FF_curr(bin == counter_bin,31:60);is_FF_curr = reshape(is_FF_curr', numel(is_FF_curr),1);
            is_FT_prev = cell2mat(count.is_FT); is_FT_prev = is_FT_prev(bin == counter_bin,1:30);is_FT_prev = reshape(is_FT_prev', numel(is_FT_prev),1);
            is_FT_curr = cell2mat(count.is_FT); is_FT_curr = is_FT_curr(bin == counter_bin,31:60);is_FT_curr = reshape(is_FT_curr', numel(is_FT_curr),1);
            is_TF_prev = cell2mat(count.is_TF); is_TF_prev = is_TF_prev(bin == counter_bin,1:30);is_TF_prev = reshape(is_TF_prev', numel(is_TF_prev),1);
            is_TF_curr = cell2mat(count.is_TF); is_TF_curr = is_TF_curr(bin == counter_bin,31:60);is_TF_curr = reshape(is_TF_curr', numel(is_TF_curr),1);

            is_T_prev = is_TT_prev | is_TF_prev;
            is_T_curr = is_TT_curr | is_TF_curr;
            is_F_prev = is_FF_prev | is_FT_prev;
            is_F_curr = is_FF_curr | is_FT_curr;

            % remove nans
            is_nan = isnan(tag_lick);
            data(is_nan) = [];
            tag_lick(is_nan) = [];
            is_T_prev(is_nan) = []; is_T_curr(is_nan) = [];
            is_F_prev(is_nan) = []; is_F_curr(is_nan) = [];
            is_TT_prev(is_nan) = []; is_TT_curr(is_nan) = [];
            is_FF_prev(is_nan) = []; is_FF_curr(is_nan) = [];
            is_FT_prev(is_nan) = []; is_FT_curr(is_nan) = [];
            is_TF_prev(is_nan) = []; is_TF_curr(is_nan) = [];

            mean_data = [nanmean(data(is_T_prev));nanmean(data(is_T_curr));...
                nanmean(data(is_F_prev));nanmean(data(is_F_curr));...
                nanmean(data(is_T_curr)-data(is_T_prev)); nanmean(data(is_F_curr)-data(is_F_prev))] ;
            sem_data = [nanstd(data(is_T_prev))/sqrt(sum(is_T_prev));nanstd(data(is_T_curr))/sqrt(sum(is_T_curr)); ...
                nanstd(data(is_F_prev))/sqrt(sum(is_F_prev));nanstd(data(is_F_curr))/sqrt(sum(is_F_curr));...
                nanstd(data(is_T_curr)-data(is_T_prev))/sqrt(sum(is_T_curr+is_T_prev)); nanstd(data(is_F_curr)-data(is_F_prev))/sqrt(sum(is_F_curr+is_F_prev))] ;

            % t test
            %             [h_T,p_T] = ttest(data(is_T_prev), data(is_T_curr));
            %             [h_F,p_F] = ttest(data(is_F_prev), data(is_F_curr));
            %             [h_diff,p_diff] = ttest2(data(is_T_curr) - data(is_T_prev), data(is_F_curr) - data(is_F_prev));
            [h_T_diff,p_T_diff,CI_T_diff,stat_T_diff] = ttest(data(is_T_curr) - data(is_T_prev));
            writetable(array2table(data(is_T_curr) - data(is_T_prev)), [path.out_path 'data' filesep path.path_data_monkey_sorted 'T_diff' ext '.csv']);

            [h_F_diff,p_F_diff,CI_F_diff,stat_F_diff] = ttest(data(is_F_curr) - data(is_F_prev));
            writetable(array2table(data(is_F_curr) - data(is_F_prev)), [path.out_path 'data' filesep path.path_data_monkey_sorted 'F_diff' ext '.csv']);

            [h_TF_diff,p_TF_diff,CI_TF_diff,stat_TF_diff] = ttest2(data(is_T_curr) - data(is_T_prev), data(is_F_curr) - data(is_F_prev));
            data_TF_diff = nan(max([sum(is_T_curr) sum(is_F_curr)]),2); data_TF_diff(1:sum(is_T_curr),1) = (data(is_T_curr) - data(is_T_prev))'; data_TF_diff(1:sum(is_F_curr),2) = (data(is_F_curr) - data(is_F_prev))'; 
            writetable(array2table(data_TF_diff), [path.out_path 'data' filesep path.path_data_monkey_sorted 'TF_diff' ext '.csv']);

            mean_T_diff = mean_data(end-1);
            SE_T_diff = sem_data(end-1);
            mean_F_diff = mean_data(end);
            SE_F_diff = sem_data(end);

            test = ['T_diff: t(' num2str(stat_T_diff.df) ') = ' num2str(stat_T_diff.tstat) ', p = '  num2str(p_T_diff) ', ' num2str(stat_T_diff.sd) ', ' num2str(mean_T_diff) '+/-' num2str(SE_T_diff) newline ...
                'F_diff: t(' num2str(stat_F_diff.df) ') = ' num2str(stat_F_diff.tstat) ', p = '  num2str(p_F_diff) ', ' num2str(stat_F_diff.sd) ', ' num2str(mean_F_diff) '+/-' num2str(SE_F_diff) newline ...
                'T-F_diff: t(' num2str(stat_TF_diff.df) ') = ' num2str(stat_TF_diff.tstat) ', p = '  num2str(p_TF_diff) ', ' num2str(stat_TF_diff.sd)];

            %%% Data for Plot (1.1)  %%%
            if contains(bin_list(counter_bin_list), 'combined') && counter_variable == 1
                for counter_tag_lick = 1 : count.num_tag_lick
                    if counter_tag_lick == 1
                        is_tag_lick = tag_lick == 1;
                    elseif counter_tag_lick == 2
                        is_tag_lick = tag_lick == 2 | tag_lick == 3;
                    elseif counter_tag_lick == 4
                        is_tag_lick = tag_lick == 4 | tag_lick == 5;
                    elseif counter_tag_lick == 6
                        is_tag_lick = tag_lick == 6 | tag_lick == 7;
                    else
                        is_tag_lick = tag_lick == 0;
                    end
                    prob_is_T_prev(1,counter_tag_lick) = nansum(is_T_prev & is_tag_lick)/length(is_tag_lick);
                    prob_is_T_curr(1,counter_tag_lick) = nansum(is_T_curr & is_tag_lick)/length(is_tag_lick);
                    prob_is_F_prev(1,counter_tag_lick) = nansum(is_F_prev & is_tag_lick)/length(is_tag_lick);
                    prob_is_F_curr(1,counter_tag_lick) = nansum(is_F_curr & is_tag_lick)/length(is_tag_lick);
                end
                inds_nan = [3 5 7];
                prob_is_T_prev(inds_nan) = nan;
                prob_is_T_curr(inds_nan) = nan;
                prob_is_F_prev(inds_nan) = nan;
                prob_is_F_curr(inds_nan) = nan;
            end
            %%% Data for Plot (1.2)  %%%
            if contains(bin_list(counter_bin_list), 'combined') && counter_variable == 1
                prob_is_TX_prev = nansum(is_T_prev)/length(is_T_prev);
                prob_is_FX_prev = nansum(is_F_prev)/length(is_F_prev);
                prob_is_TT_curr = nansum(is_TT_curr)/length(is_TT_curr);
                prob_is_FT_curr = nansum(is_FT_curr)/length(is_FT_curr);
                prob_is_TF_curr = nansum(is_TF_curr)/length(is_TF_curr);
                prob_is_FF_curr = nansum(is_FF_curr)/length(is_FF_curr);
            end

            yyaxis right
            errorbar([5 6],mean_data([5 6]), sem_data([5 6]), ['.' colors(counter_bin)], 'LineWidth', 0.75, 'MarkerSize', 1);
            plot([5 6],mean_data([5 6]), ['-' colors(counter_bin)], 'LineWidth', 0.75);
            set(gca, 'Ycolor','k')
            yyaxis left
            errorbar([1 2],mean_data([1 2]), sem_data([1 2]), ['.' colors(counter_bin)], 'LineWidth', 0.75, 'MarkerSize', 1);
            plot([1 2],mean_data([1 2]), ['-' colors(counter_bin)], 'LineWidth', 0.75);
            errorbar([3 4],mean_data([3 4]), sem_data([3 4]), ['.' colors(counter_bin)], 'LineWidth', 0.75, 'MarkerSize', 1);
            plot([3 4],mean_data([3 4]), ['-' colors(counter_bin)], 'LineWidth', 0.75);
            xline(2.5, 'k', 'LineWidth', 0.5)
            xline(4.5, 'k', 'LineWidth', 0.5)
            set(gca, 'Ycolor','k')
        end
        xlim([0.5 6.5])
        xticks([1:6])
        xticklabels(x_label)
        ylabel(variable_list(counter_variable),'Interpreter', 'none')
        ylim(y_lim)
        yticks(y_ticks)
        title_ = bin_list(counter_bin_list);
        if contains(bin_list(counter_bin_list), 'combined')
            title([title_ test ], 'Interpreter', 'none')
        else
            title(title_,'Interpreter', 'none')
        end

    end
end
ESN_Beautify_Plot(gcf, [20,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_RPE' ], 'pdf');
%close(gcf)

%%% PLOT(1.1) %%%
figure
hold on
for counter_tag_lick = 1 : count.num_tag_lick
    plot([1 2], [prob_is_T_prev(counter_tag_lick) prob_is_T_curr(counter_tag_lick)], colors(counter_tag_lick),'LineWidth', 1)
    plot([1 2], [prob_is_T_prev(counter_tag_lick) prob_is_T_curr(counter_tag_lick)], ['o' colors(counter_tag_lick)],'MarkerSize', 5)
    plot([3 4], [prob_is_F_prev(counter_tag_lick) prob_is_F_curr(counter_tag_lick)], colors(counter_tag_lick),'LineWidth', 1)
    plot([3 4], [prob_is_F_prev(counter_tag_lick) prob_is_F_curr(counter_tag_lick)],  ['o' colors(counter_tag_lick)],'MarkerSize', 5)
end
xline(2.5, 'k', 'LineWidth', 0.5)
set(gca, 'Ycolor','k')
xlim([0.5 4.5])
xticks([1:4])
xticklabels(x_label)
ylabel('probability')
ylim(y_lim)
yticks(y_ticks)

ESN_Beautify_Plot(gcf, [5,5])
saveas(gcf, [path.out_path path_ext filesep path.path_data_monkey_sorted '_RPE_type' ], 'pdf');
%close(gcf)


%%% PLOT(1.2) %%%
figure
hold on
plot([1 2], [prob_is_TX_prev prob_is_TT_curr], colors(1),'LineWidth', 1)
plot([1 2], [prob_is_TX_prev prob_is_TT_curr], ['o' colors(1)],'MarkerSize', 5)
plot([1 2], [prob_is_TX_prev prob_is_TF_curr], colors(2),'LineWidth', 1)
plot([1 2], [prob_is_TX_prev prob_is_TF_curr], ['o' colors(2)],'MarkerSize', 5)
plot([3 4], [prob_is_FX_prev prob_is_FT_curr], colors(1),'LineWidth', 1)
plot([3 4], [prob_is_FX_prev prob_is_FT_curr],  ['o' colors(1)],'MarkerSize', 5)
plot([3 4], [prob_is_FX_prev prob_is_FF_curr], colors(2),'LineWidth', 1)
plot([3 4], [prob_is_FX_prev prob_is_FF_curr],  ['o' colors(2)],'MarkerSize', 5)
xline(2.5, 'k', 'LineWidth', 0.5)
set(gca, 'Ycolor','k')
xlim([0.5 4.5])
xticks([1:4])
xticklabels(x_label)
ylabel('probability')
ylim(y_lim)
yticks(y_ticks)

ESN_Beautify_Plot(gcf, [5,5])
saveas(gcf, [path.out_path path_ext filesep path.path_data_monkey_sorted '_RPE_success' ], 'pdf');
%close(gcf)

%% FIGURE - lick prob
clearvars -except params lick sac path count colors count_session ext path_ext plot_session plot_stat
path_ext =  ['HARVEST' filesep ];

% compute % lick
for counter_session = 1 : count.num_session
    tag_lick = lick.tag_lick{counter_session,1};tag_lick = reshape(tag_lick', numel(tag_lick),1);tag_lick(isnan(tag_lick)) = [];
    prob_inner_tube(counter_session,1) = sum(tag_lick == 2 | tag_lick == 3)/sum(tag_lick == 2 | tag_lick == 3 | tag_lick == 5 | tag_lick == 7) * 100;
end
mean_prob_inner_tube = nanmean(prob_inner_tube);
sem_prob_inner_tube = nanstd(prob_inner_tube)/sqrt(size(prob_inner_tube,1));

print_ = 0;
if print_ == 1
    disp([num2str(mean_prob_inner_tube) ' +/- ' num2str(sem_prob_inner_tube)])
end

%%% PLOT (1) %%%
bin_list = ["interval_bin_", "rew_available_bin_", "dist_tube_bin_", "num_trial_bin_"];
num_type = 2;
subplot_even = [2:2:length(bin_list)*num_type];
subplot_odd = [1:2:length(bin_list)*num_type];

figure
for counter_bin_list = 1 : length(bin_list)

    %rew driven vs groom licks
    subplot(length(bin_list),num_type,subplot_odd(counter_bin_list))
    prob_tag_lick_rew = nan(count.num_session,5);
    mean_prob_tag_lick_rew = [];
    se_prob_tag_lick_rew = [];
    prob_tag_lick_groom = nan(count.num_session,5);
    mean_prob_tag_lick_groom = [];
    se_prob_tag_lick_groom = [];
    %     bin = cell2mat(count.(bin_list(counter_bin_list)));
    for counter_session = 1 : count.num_session
        bin = count.(bin_list(counter_bin_list)){counter_session,1};
        for counter_bin = 1 : max(bin)
            tag_lick = lick.tag_lick{counter_session,1};tag_lick = tag_lick(bin == counter_bin,:);tag_lick = reshape(tag_lick', numel(tag_lick),1);tag_lick(isnan(tag_lick)) = [];
            if isempty(tag_lick)
                prob_tag_lick_rew(counter_session,counter_bin) = nan;
                prob_tag_lick_groom(counter_session,counter_bin) = nan;
            else
                prob_tag_lick_rew(counter_session,counter_bin) = nansum(tag_lick ~=1)/length(tag_lick) ;
                prob_tag_lick_groom(counter_session,counter_bin) = nansum(tag_lick == 1)/length(tag_lick);
            end
        end

    end
    mean_prob_tag_lick_rew = nanmean(prob_tag_lick_rew);
    se_prob_tag_lick_rew = nanstd(prob_tag_lick_rew)/sqrt(count.num_session);
    mean_prob_tag_lick_groom = nanmean(prob_tag_lick_groom);
    se_prob_tag_lick_groom = nanstd(prob_tag_lick_groom)/sqrt(count.num_session);
    hold on
    errorbar(1:length(mean_prob_tag_lick_rew), mean_prob_tag_lick_rew, se_prob_tag_lick_rew, ['.' colors(1)], 'LineWidth', 2, 'MarkerSize', 10);
    plot(1:length(mean_prob_tag_lick_rew), mean_prob_tag_lick_rew, ['-' colors(1)], 'LineWidth',1)
    ylabel('prob. rew lick')
    yyaxis right;
    errorbar(1:length(mean_prob_tag_lick_groom), mean_prob_tag_lick_groom, se_prob_tag_lick_groom, ['.' colors(2)], 'LineWidth', 2, 'MarkerSize', 10);
    plot(1:length(mean_prob_tag_lick_groom), mean_prob_tag_lick_groom, ['-' colors(2)], 'LineWidth',1)
    ylabel('prob. groom lick')
    set(gca, 'Ycolor',colors(2))
    xlabel(bin_list(counter_bin_list),'Interpreter', 'none')
    xlim([0.5 length(mean_prob_tag_lick_rew)+0.5])
    xticks([1:length(mean_prob_tag_lick_rew)])
    title(['% inner tube from rew seeking: ' num2str(mean_prob_inner_tube) ' +/- ' num2str(sem_prob_inner_tube)])

    % rew driven: T vs F
    subplot(length(bin_list),num_type,subplot_even(counter_bin_list))
    prob_tag_lick_T = nan(count.num_session,5);
    mean_prob_tag_lick_T = [];
    se_prob_tag_lick_T = [];
    prob_tag_lick_F = nan(count.num_session,5);
    mean_prob_tag_lick_F = [];
    se_prob_tag_lick_F = [];
    % bin = cell2mat(count.(bin_list(counter_bin_list)));
    for counter_session = 1 : count.num_session
        bin = count.(bin_list(counter_bin_list)){counter_session,1};
        for counter_bin = 1 : max(bin)
            tag_lick = lick.tag_lick{counter_session,1}; tag_lick = tag_lick(bin == counter_bin,:);tag_lick = reshape(tag_lick', numel(tag_lick),1);tag_lick(isnan(tag_lick)|tag_lick == 1|tag_lick == 4|tag_lick == 6) = [];
            if isempty(tag_lick)
                prob_tag_lick_T(counter_session,counter_bin) = nan;
                prob_tag_lick_F(counter_session,counter_bin) = nan;
            else
                prob_tag_lick_T(counter_session,counter_bin) = nansum(~mod(tag_lick,2) & tag_lick~= 3)/length(tag_lick) ;
                prob_tag_lick_F(counter_session,counter_bin) = nansum(mod(tag_lick,2) | tag_lick == 3)/length(tag_lick) ;
            end
        end
    end
    mean_prob_tag_lick_T = nanmean(prob_tag_lick_T);
    se_prob_tag_lick_T = nanstd(prob_tag_lick_T)/sqrt(count.num_session);
    mean_prob_tag_lick_F = nanmean(prob_tag_lick_F);
    se_prob_tag_lick_F = nanstd(prob_tag_lick_F)/sqrt(count.num_session);

    hold on
    errorbar(1:length(mean_prob_tag_lick_T), mean_prob_tag_lick_T, se_prob_tag_lick_T, ['.' colors(1)], 'LineWidth', 2, 'MarkerSize', 10);
    plot(1:length(mean_prob_tag_lick_T), mean_prob_tag_lick_T, ['-' colors(1)], 'LineWidth',1)
    ylabel('prob. T lick')
    yyaxis right;
    errorbar(1:length(mean_prob_tag_lick_F), mean_prob_tag_lick_F, se_prob_tag_lick_F, ['.' colors(2)], 'LineWidth', 2, 'MarkerSize', 10);
    plot(1:length(mean_prob_tag_lick_F), mean_prob_tag_lick_F, ['-' colors(2)], 'LineWidth',1)
    ylabel('prob. F lick')
    set(gca, 'Ycolor',colors(2))
    xlabel(bin_list(counter_bin_list),'Interpreter', 'none')
    xlim([0.5 length(mean_prob_tag_lick_T)+0.5])
    xticks([1:length(mean_prob_tag_lick_T)])
    title('')
end
ESN_Beautify_Plot(gcf, [10,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_prob_Xbin' ], 'pdf');
close(gcf)

%%% PLOT (2) %%%
bin_list = ["interval_bin_", "rew_available_bin_", "dist_tube_bin_", "num_trial_bin_"];
num_type = 2;
subplot_even = [2:2:length(bin_list)*num_type];
subplot_odd = [1:2:length(bin_list)*num_type];

figure
for counter_bin_list = 1 : length(bin_list)

    %rew driven vs groom licks
    subplot(length(bin_list),num_type,subplot_odd(counter_bin_list))
    prob_tag_lick_rew = [];
    se_prob_tag_lick_rew = [];
    prob_tag_lick_groom = [];
    se_prob_tag_lick_groom = [];
    hold on
    bin = cell2mat(count.(bin_list(counter_bin_list)));
    for counter_bin = 1 : max(bin)
        tag_lick = cell2mat(lick.tag_lick); tag_lick = tag_lick(bin == counter_bin,:);
        prob_tag_lick_rew(counter_bin,:) = nansum(tag_lick ~=1)/size(tag_lick,1);
        se_prob_tag_lick_rew(counter_bin,:) = nanstd(tag_lick ~=1)/size(tag_lick,1);
        prob_tag_lick_groom(counter_bin,:) = nansum(tag_lick == 1)/size(tag_lick,1);
        se_prob_tag_lick_groom(counter_bin,:) = nanstd(tag_lick == 1)/size(tag_lick,1);

        yyaxis left;
        errorbar(1:length(prob_tag_lick_rew), prob_tag_lick_rew(counter_bin,:), se_prob_tag_lick_rew(counter_bin,:), [ '.' colors(counter_bin)], 'LineWidth', 2, 'MarkerSize', 10);
        %shade(1:length(prob_tag_lick_rew), prob_tag_lick_rew(counter_bin,:) + se_prob_tag_lick_rew(counter_bin,:), [ '-' colors(counter_bin)], 1:length(prob_tag_lick_rew), prob_tag_lick_rew(counter_bin,:) - se_prob_tag_lick_rew(counter_bin,:),[ '-' colors(counter_bin)], 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',[colors(counter_bin)], 'FillColor', colors(counter_bin));
        plot(1:length(prob_tag_lick_rew), prob_tag_lick_rew(counter_bin,:), [ '-' colors(counter_bin)], 'LineWidth',1)
        yyaxis right;
        errorbar(1:length(prob_tag_lick_groom), prob_tag_lick_groom(counter_bin,:), se_prob_tag_lick_groom(counter_bin,:), [ '.' colors(counter_bin)], 'LineWidth', 2, 'MarkerSize', 10);
        %shade(1:length(prob_tag_lick_groom), prob_tag_lick_groom(counter_bin,:) + se_prob_tag_lick_groom(counter_bin,:), [ '--' colors(counter_bin)],1:length(prob_tag_lick_groom), prob_tag_lick_groom(counter_bin,:) - se_prob_tag_lick_groom(counter_bin,:),[ '--' colors(counter_bin)], 'LineWidth', 0.1, 'FillType',[1 2;2 1],'Color',[colors(counter_bin)], 'FillColor', colors(counter_bin));
        plot(1:length(prob_tag_lick_groom), prob_tag_lick_groom(counter_bin,:), ['--' colors(counter_bin)], 'LineWidth',1)
    end
    yyaxis left;
    ylabel('prob. rew lick')
    set(gca, 'Ycolor',colors(1))
    ylim([0.85 1])

    yyaxis right;
    ylabel('prob. groom lick')
    set(gca, 'Ycolor',colors(1))
    ylim([0 0.15])
    xlabel('num lick in harvest')
    xlim([0.5 length(prob_tag_lick_rew)+0.5])
    xticks([1 10 20 30])
    title_ = bin_list(counter_bin_list);
    title(['All licks: -reward seeking, --grooming | ' title_],'Interpreter', 'none')

    % rew driven: T vs F
    subplot(length(bin_list),num_type,subplot_even(counter_bin_list))
    prob_tag_lick_T = [];
    se_prob_tag_lick_T = [];
    prob_tag_lick_F = [];
    se_prob_tag_lick_F = [];
    hold on
    bin = cell2mat(count.(bin_list(counter_bin_list)));
    for counter_bin = 1 : max(bin)
        tag_lick = cell2mat(lick.tag_lick); tag_lick = tag_lick(bin == counter_bin,:);tag_lick(isnan(tag_lick))=1; % make nans 1, they will be disregarded
        prob_tag_lick_T(counter_bin,:) = nansum(~mod(tag_lick,2) & tag_lick ~= 1)/size(tag_lick,1);
        se_prob_tag_lick_T(counter_bin,:) = nanstd(~mod(tag_lick,2) & tag_lick ~= 1)/size(tag_lick,1);
        prob_tag_lick_F(counter_bin,:) = nansum(mod(tag_lick,2) & tag_lick ~= 1)/size(tag_lick,1);
        se_prob_tag_lick_F(counter_bin,:) = nanstd(mod(tag_lick,2) & tag_lick ~= 1)/size(tag_lick,1);

        yyaxis left;
        errorbar(1:length(prob_tag_lick_T), prob_tag_lick_T(counter_bin,:), se_prob_tag_lick_T(counter_bin,:), [ '.' colors(counter_bin)], 'LineWidth', 2, 'MarkerSize', 10);
        plot(1:length(prob_tag_lick_T), prob_tag_lick_T(counter_bin,:), [ '-' colors(counter_bin)], 'LineWidth',1)
        yyaxis right;
        errorbar(1:length(prob_tag_lick_F), prob_tag_lick_F(counter_bin,:), se_prob_tag_lick_F(counter_bin,:), [ '.' colors(counter_bin)], 'LineWidth', 2, 'MarkerSize', 10);
        plot(1:length(prob_tag_lick_F), prob_tag_lick_F(counter_bin,:), ['--' colors(counter_bin)], 'LineWidth',1)
    end
    yyaxis left;
    ylabel('prob. T lick')
    set(gca, 'Ycolor',colors(1))
    %ylim([0.85 1])

    yyaxis right;
    ylabel('prob. F lick')
    set(gca, 'Ycolor',colors(1))
    %ylim([0 0.15])
    xlabel('num lick in harvest')
    xlim([0.5 length(prob_tag_lick_T)+0.5])
    xticks([1 10 20 30])
    title_ = bin_list(counter_bin_list);
    title(['Reward seeking licks -T, --F | ' title_],'Interpreter', 'none')
end
ESN_Beautify_Plot(gcf, [10,10])
saveas(gcf, [path.out_path path_ext path.path_data_monkey_sorted '_lick_prob_Xlick' ], 'pdf');
close(gcf)



%% Figure tongue kinematics
% clearvars -except params lick sac path count
%
%
% %%% concatenate tongue_vm %%%
% for counter_session = 1 : count.num_session
%     tongue_vm{counter_session,1} = cell2mat(lick.tongue_vm{counter_session,1});
%     tongue_dm{counter_session,1} = cell2mat(lick.tongue_dm{counter_session,1});
%     tongue_ang{counter_session,1} = cell2mat(lick.tongue_ang{counter_session,1});
%     tongue_tip_px{counter_session,1} = cell2mat(lick.tongue_tip_px{counter_session,1});
%     tongue_tip_py{counter_session,1} = cell2mat(lick.tongue_tip_py{counter_session,1});
% end
% tongue_vm = cell2mat(tongue_vm);
% is_nan = isnan(tongue_vm(:,1));
% tongue_vm(is_nan,:) = [];
% tongue_dm = cell2mat(tongue_dm);tongue_dm(is_nan,:) = [];
% tongue_ang = cell2mat(tongue_ang);tongue_ang(is_nan,:) = [];
% tongue_tip_px = cell2mat(tongue_tip_px);tongue_tip_px(is_nan,:) = [];
% tongue_tip_py = cell2mat(tongue_tip_py);tongue_tip_py(is_nan,:) = [];
% dist_err_tongue_rew  = cell2mat(lick.dist_err_tongue_rew); dist_err_tongue_rew = reshape(dist_err_tongue_rew', numel(dist_err_tongue_rew),1);dist_err_tongue_rew(is_nan,:) = [];
% tongue_tip_px_dmax = cell2mat(lick.tongue_tip_px_dmax); tongue_tip_px_dmax = reshape(tongue_tip_px_dmax', numel(tongue_tip_px_dmax),1);tongue_tip_px_dmax(is_nan,:) = [];
% tongue_tip_py_dmax = cell2mat(lick.tongue_tip_py_dmax); tongue_tip_py_dmax = reshape(tongue_tip_py_dmax', numel(tongue_tip_py_dmax),1);tongue_tip_py_dmax(is_nan,:) = [];
% tongue_mid_px_dmax = cell2mat(lick.tongue_mid_px_dmax); tongue_mid_px_dmax = reshape(tongue_mid_px_dmax', numel(tongue_mid_px_dmax),1);tongue_mid_px_dmax(is_nan,:) = [];
% tongue_mid_py_dmax = cell2mat(lick.tongue_mid_py_dmax); tongue_mid_py_dmax = reshape(tongue_mid_py_dmax', numel(tongue_mid_py_dmax),1);tongue_mid_py_dmax(is_nan,:) = [];
% tongue_r_px_dmax = cell2mat(lick.tongue_r_px_dmax); tongue_r_px_dmax = reshape(tongue_r_px_dmax', numel(tongue_r_px_dmax),1);tongue_r_px_dmax(is_nan,:) = [];
% tongue_r_py_dmax = cell2mat(lick.tongue_r_py_dmax); tongue_r_py_dmax = reshape(tongue_r_py_dmax', numel(tongue_r_py_dmax),1);tongue_r_py_dmax(is_nan,:) = [];
% tongue_l_px_dmax = cell2mat(lick.tongue_l_px_dmax); tongue_l_px_dmax = reshape(tongue_l_px_dmax', numel(tongue_l_px_dmax),1);tongue_l_px_dmax(is_nan,:) = [];
% tongue_l_py_dmax = cell2mat(lick.tongue_l_py_dmax); tongue_l_py_dmax = reshape(tongue_l_py_dmax', numel(tongue_l_py_dmax),1);tongue_l_py_dmax(is_nan,:) = [];
% rew_r_px_dmax = cell2mat(lick.rew_r_px_dmax); rew_r_px_dmax = reshape(rew_r_px_dmax', numel(rew_r_px_dmax),1);rew_r_px_dmax(is_nan,:) = [];
% rew_r_py_dmax = cell2mat(lick.rew_r_py_dmax); rew_r_py_dmax = reshape(rew_r_py_dmax', numel(rew_r_py_dmax),1);rew_r_py_dmax(is_nan,:) = [];
% rew_l_px_dmax = cell2mat(lick.rew_l_px_dmax); rew_l_px_dmax = reshape(rew_l_px_dmax', numel(rew_l_px_dmax),1);rew_l_px_dmax(is_nan,:) = [];
% rew_l_py_dmax = cell2mat(lick.rew_l_py_dmax); rew_l_py_dmax = reshape(rew_l_py_dmax', numel(rew_l_py_dmax),1);rew_l_py_dmax(is_nan,:) = [];
% nose_r_px_dmax = cell2mat(lick.nose_r_px_dmax); nose_r_px_dmax = reshape(nose_r_px_dmax', numel(nose_r_px_dmax),1);nose_r_px_dmax(is_nan,:) = [];
% nose_r_py_dmax = cell2mat(lick.nose_r_py_dmax); nose_r_py_dmax = reshape(nose_r_py_dmax', numel(nose_r_py_dmax),1);nose_r_py_dmax(is_nan,:) = [];
% nose_l_px_dmax = cell2mat(lick.nose_l_px_dmax); nose_l_px_dmax = reshape(nose_l_px_dmax', numel(nose_l_px_dmax),1);nose_l_px_dmax(is_nan,:) = [];
% nose_l_py_dmax = cell2mat(lick.nose_l_py_dmax); nose_l_py_dmax = reshape(nose_l_py_dmax', numel(nose_l_py_dmax),1);nose_l_py_dmax(is_nan,:) = [];
% rtube_r_px_dmax = cell2mat(lick.rtube_r_px_dmax); rtube_r_px_dmax = reshape(rtube_r_px_dmax', numel(rtube_r_px_dmax),1);rtube_r_px_dmax(is_nan,:) = [];
% rtube_r_py_dmax = cell2mat(lick.rtube_r_py_dmax); rtube_r_py_dmax = reshape(rtube_r_py_dmax', numel(rtube_r_py_dmax),1);rtube_r_py_dmax(is_nan,:) = [];
% rtube_l_px_dmax = cell2mat(lick.rtube_l_px_dmax); rtube_l_px_dmax = reshape(rtube_l_px_dmax', numel(rtube_l_px_dmax),1);rtube_l_px_dmax(is_nan,:) = [];
% rtube_l_py_dmax = cell2mat(lick.rtube_l_py_dmax); rtube_l_py_dmax = reshape(rtube_l_py_dmax', numel(rtube_l_py_dmax),1);rtube_l_py_dmax(is_nan,:) = [];
% ltube_r_px_dmax = cell2mat(lick.ltube_r_px_dmax); ltube_r_px_dmax = reshape(ltube_r_px_dmax', numel(ltube_r_px_dmax),1);ltube_r_px_dmax(is_nan,:) = [];
% ltube_r_py_dmax = cell2mat(lick.ltube_r_py_dmax); ltube_r_py_dmax = reshape(ltube_r_py_dmax', numel(ltube_r_py_dmax),1);ltube_r_py_dmax(is_nan,:) = [];
% ltube_l_px_dmax = cell2mat(lick.ltube_l_px_dmax); ltube_l_px_dmax = reshape(ltube_l_px_dmax', numel(ltube_l_px_dmax),1);ltube_l_px_dmax(is_nan,:) = [];
% ltube_l_py_dmax = cell2mat(lick.ltube_l_py_dmax); ltube_l_py_dmax = reshape(ltube_l_py_dmax', numel(ltube_l_py_dmax),1);ltube_l_py_dmax(is_nan,:) = [];
% tag_lick = cell2mat(lick.tag_lick); tag_lick = reshape(tag_lick', numel(tag_lick),1);tag_lick(is_nan,:) = [];
%
%
% %%% cluster protrusion phase %%%
% % for counter_lick = 1 : size(tongue_vm,1)
% %     tongue_vm( counter_lick,tongue_vm(counter_lick,:) < 0) = 0;
% % end
%
% % tongue_vm_fail = tongue_vm(tag_lick==5,:);
% %tongue_vm = tongue_vm_fail;
%
% peaks = nan(size(tongue_vm,1), size(tongue_vm,2));
% is_double_peak = logical(zeros(size(tongue_vm,1),1));
% for counter_lick = 1 : size(tongue_vm,1)
%     [peaks_,~,~,~] = findpeaks(tongue_vm(counter_lick,:), 'MinPeakProminence',max(tongue_vm(counter_lick,:))*0.80 , 'MinPeakDistance',50,'Annotate', 'extent');
%     peaks(counter_lick,1:length(peaks_)) = peaks_;
%     if length(peaks_) > 1
%         is_double_peak(counter_lick) = 1;
%     end
% end
% % VIEWER
% %%% viewer %%%
% %n_lick = find(tag_lick == 6);
% n_lick = 7;
%
% % figure
% % histogram(tag_lick(is_double_peak),'FaceColor', 'k')
% % xticklabels({"Grm", "IT_T", "IT_F", "OE_T", "OE_F", "UT_T", "UT_F"})
% % xticks(1:length(xticklabels))
% % xlabel('lick type')
% % ylabel('count')
% % title('lick types with double peak pro. velocity profile; 80% prom.')
% % ESN_Beautify_Plot(gcf,[8,8])
% % saveas(gcf, [path.out_path 'LICK' filesep path.path_data_monkey_sorted '_tongue_vm_hist' ], 'pdf');
%
% for counter_lick = 1 : length(n_lick)
%     figure
%     subplot(3,2,1)
%     plot(tongue_dm(n_lick(counter_lick),:),'k','LineWidth', 2);
%     xlim([-50 500])
%     ylabel('displacement (mm)')
%     title(['Lick: ' num2str(n_lick(counter_lick)) ' | Tag: ' num2str(tag_lick(n_lick(counter_lick))) ' | Error: ' num2str(dist_err_tongue_rew(n_lick(counter_lick),:)) 'mm'])
%     subplot(3,2,3)
%     hold on
%     plot(tongue_vm(n_lick(counter_lick),:),'k','LineWidth', 2);
%     xlabel('time (ms)')
%     ylabel('velocity (mm/s)')
%     xlim([-50 500])
%
%     subplot(3,2,5)
%     hold on
%     plot(tongue_ang(n_lick(counter_lick),:),'k','LineWidth', 2);
%     xlabel('time (ms)')
%     ylabel('angle (deg)')
%     xlim([-50 500])
%
%
%     subplot(3,2,[2 4 6])
%     hold on
%     line([rtube_r_px_dmax(n_lick(counter_lick),:) rtube_r_px_dmax(n_lick(counter_lick),:)] ,[rtube_r_py_dmax(n_lick(counter_lick),:) 20], 'Color', 'black','LineWidth', 2);
%     line([rtube_l_px_dmax(n_lick(counter_lick),:) rtube_l_px_dmax(n_lick(counter_lick),:)] ,[rtube_l_py_dmax(n_lick(counter_lick),:) 20], 'Color', 'black','LineWidth', 2);
%     line([rtube_r_px_dmax(n_lick(counter_lick),:) rtube_l_px_dmax(n_lick(counter_lick),:)] , [rtube_r_py_dmax(n_lick(counter_lick),:) rtube_l_py_dmax(n_lick(counter_lick),:)], 'Color', 'black','LineWidth', 2);
%     line([ltube_r_px_dmax(n_lick(counter_lick),:) ltube_r_px_dmax(n_lick(counter_lick),:)] ,[-20 ltube_r_py_dmax(n_lick(counter_lick),:)], 'Color', 'black','LineWidth', 2);
%     line([ltube_l_px_dmax(n_lick(counter_lick),:) ltube_l_px_dmax(n_lick(counter_lick),:)] ,[-20 ltube_l_py_dmax(n_lick(counter_lick),:)], 'Color', 'black','LineWidth', 2);
%     line([ltube_r_px_dmax(n_lick(counter_lick),:) ltube_l_px_dmax(n_lick(counter_lick),:)], [ltube_r_py_dmax(n_lick(counter_lick),:) ltube_l_py_dmax(n_lick(counter_lick),:)], 'Color', 'black','LineWidth', 2);
%     plot(tongue_tip_px(n_lick(counter_lick),:),tongue_tip_py(n_lick(counter_lick),:),'k','LineWidth', 2);
%     plot(tongue_tip_px_dmax(n_lick(counter_lick),:),tongue_tip_py_dmax(n_lick(counter_lick),:),'ok');
%     plot(tongue_mid_px_dmax(n_lick(counter_lick),:),tongue_mid_py_dmax(n_lick(counter_lick),:),'oc');
%     plot(tongue_l_px_dmax(n_lick(counter_lick),:),tongue_l_py_dmax(n_lick(counter_lick),:),'ob');
%     plot(tongue_r_px_dmax(n_lick(counter_lick),:),tongue_r_py_dmax(n_lick(counter_lick),:),'or');
%     plot(rew_r_px_dmax(n_lick(counter_lick),:),rew_r_py_dmax(n_lick(counter_lick),:),'.g','MarkerSize', 15);
%     plot(rew_l_px_dmax(n_lick(counter_lick),:),rew_l_py_dmax(n_lick(counter_lick),:),'.g','MarkerSize', 15);
%     plot(nose_r_px_dmax(n_lick(counter_lick),:),nose_r_py_dmax(n_lick(counter_lick),:),'ok');
%     plot(nose_l_px_dmax(n_lick(counter_lick),:),nose_l_py_dmax(n_lick(counter_lick),:),'ok');
%     plot(rtube_r_px_dmax(n_lick(counter_lick),:),rtube_r_py_dmax(n_lick(counter_lick),:),'.y','MarkerSize', 15);
%     plot(rtube_l_px_dmax(n_lick(counter_lick),:),rtube_l_py_dmax(n_lick(counter_lick),:),'.y','MarkerSize', 15);
%     plot(ltube_r_px_dmax(n_lick(counter_lick),:),ltube_r_py_dmax(n_lick(counter_lick),:),'.y','MarkerSize', 15);
%     plot(ltube_l_px_dmax(n_lick(counter_lick),:),ltube_l_py_dmax(n_lick(counter_lick),:),'.y','MarkerSize', 15);
%     xlabel('x (mm)')
%     ylabel('y (mm)')
%     xlim([-10, 20]);
%     ylim([-20, 20]);
%     set(gca, 'YDir','reverse')
%     ESN_Beautify_Plot(gcf, [20,10])
% %     saveas(gcf, [path.out_path 'LICK' filesep path.path_data_monkey_sorted '_tongue_vm_' num2str(n_lick(counter_lick)) '_' num2str(tag_lick(n_lick(counter_lick)))], 'pdf');
%     % close(gcf)
% end

%% DR and CLUSTER
% [~, pca_mat, ~] = pca(tongue_vm);
% [reduction, umap, clusterIdentifiers, extras]=run_umap(tongue_vm);
%
% idx_pca = kmeans([pca_mat(:,1) pca_mat(:,2)],3);
% idx_umap = kmeans([reduction(:,1) reduction(:,2)],3);
%
% figure
% subplot(2,2,1)
% hold on
% plot(pca_mat(idx_pca == 1, 1), pca_mat(idx_pca == 1, 2), '.r')
% plot(pca_mat(idx_pca == 2, 1), pca_mat(idx_pca == 2, 2), '.b')
% plot(pca_mat(idx_pca == 3, 1), pca_mat(idx_pca == 3, 2), '.g')
% xlabel('PC1')
% ylabel('PC2')
%
% subplot(2,2,2)
% hold on
% plot(nanmean(tongue_vm(idx_pca == 1, :)), 'r')
% plot(nanmean(tongue_vm(idx_pca == 2, :)), 'b')
% plot(nanmean(tongue_vm(idx_pca == 3, :)), 'g')
% xlabel('time (ms)')
% ylabel('velocity pro. (mm/s)')
%
% subplot(2,2,3)
% plot(reduction(idx_umap == 1, 1), reduction(idx_umap == 1, 2), '.r')
% plot(reduction(idx_umap == 2, 1), reduction(idx_umap == 2, 2), '.b')
% plot(reduction(idx_umap == 3, 1), reduction(idx_umap == 3, 2), '.g')
% xlabel('UMAP1')
% ylabel('UMAP2')
%
% subplot(2,2,4)
% hold on
% plot(nanmean(tongue_vm(idx_umap == 1, :)), 'r')
% plot(nanmean(tongue_vm(idx_umap == 2, :)), 'b')
% plot(nanmean(tongue_vm(idx_umap == 3, :)), 'g')
% xlabel('time (ms)')
% ylabel('velocity pro. (mm/s)')
% ESN_Beautify_Plot(gcf, [20,10])
% saveas(gcf, [path.out_path 'LICK' filesep path.path_data_monkey_sorted '_cluster_tongue_vm' ], 'pdf');
% % close(gcf)
end
