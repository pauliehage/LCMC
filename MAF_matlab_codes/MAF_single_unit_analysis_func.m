%% MASTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function MAF_single_unit_analysis
function MAF_single_unit_analysis_func
%% Steps to prepare the combined eye & neuro data
clc; clear; close all;
tic;

path_data_monkey_sorted = 'data_125d';

session_list = {'2021-11-17'};

if strcmp(session_list{1}, 'all')
    session_list = dir([path_data_monkey_sorted filesep '20*' filesep '20*']);
    session_list([session_list.isdir] == 0) = [];
    session_list = {session_list.name};
end

%% params
params.sac.ang_step  = 45;
params.sac.length_trace = 500;
params.sac.inds_span        = ((-(params.sac.length_trace/2)+1) : 1 : (params.sac.length_trace/2))';
params.sac.ang_step         = 45;
params.sac.ang_edges        = (0 - (params.sac.ang_step/2)) : params.sac.ang_step : (360 + (params.sac.ang_step/2));
params.sac.ang_values       = (0) : params.sac.ang_step : (360 - params.sac.ang_step);
params.sac.tags_CS_ang_avg  = [1,4,6];
params.sac.tag_name_list    = { ...
    'prim_success', ... % tag 1
    'prim_attempt', ... % tag 2
    'prim_fail', ... % tag 3
    'corr_success', ... % tag 4
    'corr_fail', ... % tag 5
    'back_center_success', ... % tag 6
    'back_center_prim', ... % tag 7
    'back_center_irrelev', ... % tag 8
    'target_irrelev', ... % tag 9
    'other_irrelev', ... % tag 10
    'prim_no_corr',... % tag 11; prim. sac. that is not followed by corr. sac.
    'db_corr_success',... % tag 12; sac. that follows first corr. sac., back to 2nd jumped cue
    'corr_no_db_corr',... % tag 13; corr. sac. that is not followed by another corr. sac.
    'other_irrelev_visual',... % tag 14; like tag 10, but visual ang. based on the pursuit target present on the screen
    'back_center_irrelev_visual',... % tag 15; like tag 10, but visual start position based on offset of previous saccade
    ... % Add more tags here, do not reorder or change the tags defined above.
    };
params.sac.tags_CS_bundle_app = [1,4,6];
params.sac.tags_SS_bundle_app = [1,6,7];
params.sac.align_CS_bundle_app = 'visual';
params.sac.align_SS_bundle_app = 'onset';

params.lick.length_trace    = 500;
params.lick.inds_span       = ((-(params.lick.length_trace/2)+1) : 1 : (params.lick.length_trace/2))';
params.lick.ang_step        = 45;
params.lick.ang_edges       = -90 - params.lick.ang_step /2 : params.lick.ang_step  : 90 + params.lick.ang_step /2 ;
params.lick.ang_values      = -90 : params.lick.ang_step : 90;
params.lick.amp_step    = 5;
params.lick.amp_edges    = 0 : params.lick.amp_step : 25;
params.lick.vel_step    = 150;
params.lick.vel_edges    = 0 : params.lick.vel_step : 750;
params.lick.dur_step = 100;
params.lick.dur_edges = 0 : params.lick.dur_step : 500;
params.lick.tags_CS_ang_avg  = [1:9];
params.lick.tag_name_list = {  ...
    'groom', ... % tag 1
    'inner_tube_success', ... % tag 2
    'inner_tube_fail', ... % tag 3
    'outer_tube_success' ..., % tag 4
    'outer_tube_fail', ... % tag 5
    'bout_start', ... % tag 6
    'bout_end', ... % tag 7
    'harvest_start',  ...% tag 8
    'harvest_end', ... % tag 9
    };
params.lick.tag_bout_name_list = {  ...
    'bout_start_r', ... % tag 1
    'bout_start_l', ... % tag 2
    'bout_start_g', ... % tag 3
    'bout_end_r' ..., % tag 4
    'bout_end_l', ... % tag 5
    'bout_end_g', ... % tag 6
    };
params.lick.tags_CS_bundle_app = 1:9;
params.lick.tags_SS_bundle_app = 1:9;
params.lick.align_CS_bundle_app = 'onset';
params.lick.align_SS_bundle_app = 'onset';

%% funcs
funcs.file_sorter           = @PGH_reorganize_files;
funcs.spike_sorter          = @KS;
funcs.monkey_behavior_sac   = @MAF_monkey_behavior_sac;
funcs.Sac_Sorter            = @JSP_Sac_Sorter;
funcs.Sac_Recal             = @JSP_recalibrate_all_sac;
funcs.add_ephys_sac_sorter  = @MAF_add_ephys_sac_sorter;
funcs.sac_cs_on_analysis    = @JSP_CS_on_analysis_fr;
funcs.buildSacData          = @ESN_buildSacData;

funcs.monkey_behavior_lick  = @PGH_monkey_behavior_lick;
funcs.add_ephys_lick_sorter = @PGH_add_ephys_lick_sorter;
funcs.lick_cs_on_analysis   = @PGH_CS_on_analysis;
funcs.buildLickData         = @PGH_buildLickData;
funcs.buildBoutData         = @PGH_buildBoutData;

funcs.unit_summary_func     = @MAF_unit_summary_func;


%% Reorganize Files
% reorganize_files(path_data_monkey_sorted, session_list, params, funcs);

%% Spike Sorting
% sorter(path_data_monkey_sorted, session_list, params, funcs); % this func sorts

%% Behavior Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sac_analysis(path_data_monkey_sorted, session_list, params, funcs); % this func converts and runs sac_sorter on each recording .fhd files
% recalibrate(path_data_monkey_sorted, session_list, params, funcs); % this func recalibrates the tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tongue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLC(path_data_monkey_sorted, session_list, params, funcs); % this func analyzes videos and produces csv files
lick_analysis(path_data_monkey_sorted, session_list, params,funcs); % this func runs lick_sorter and produces LICKS_ALL data from .csv files

%% MetaData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% meta data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate_session_meta_data(path_data_monkey_sorted, session_list, params, funcs); % this func creates metadata

%% Unit Extraction
unit_extraction(path_data_monkey_sorted, session_list, params, funcs); % this func extracts units from sorter results

%% Alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sac_alignment(path_data_monkey_sorted, session_list, params, funcs); % this func aligns eye with ephys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tongue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lick_alignment(path_data_monkey_sorted, session_list, params, funcs); % this func aligns tongue (video) with ephys

%% Neuro Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sac_neuro_modulation(path_data_monkey_sorted, session_list, params, funcs); % this func generates combined ephys sac files in _sac.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tongue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lick_neuro_modulation(path_data_monkey_sorted, session_list, params, funcs); % this func generates combined ephys lick files in _lick.mat

%% Bundling
generate_rec_units_summary(path_data_monkey_sorted, session_list, params, funcs); % this function generates rec units summary and metadata
MAF_bundling_app(path_data_monkey_sorted, session_list{1}, params, funcs); % opens bundling app for the specified session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eye %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine_bundle_files(path_data_monkey_sorted, session_list, params, funcs);
%% Summary App
% extract_cell_interactions(path_data_monkey_sorted, session_list);
% MAF_neuromodulation_app(path_data_monkey_sorted, session_list{1},params,funcs);

toc

end

%% function reorganize files
function reorganize_files(path_data_monkey_sorted, sess_list, params)

if ~strcmp(path_data_monkey_sorted(end), filesep)
    path_data_monkey_sorted = [path_data_monkey_sorted filesep];
end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' rec ', current_rec ' \n']);
        rec_path = [sess_path current_rec filesep];
        func.file_sorter(rec_path, params, func);
    end

end
fprintf('### ALL DONE. ###\n')
end

%% function sorter
function sorter(path_data_monkey_sorted, sess_list, params, funcs)

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    geom = input('Checker (0), Linear(1)?\n');

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' \n']);
        path_to_raw = [sess_path current_rec filesep 'raw_data'];
        if exist(path_to_raw,'dir')
            try
                spike_sorter(path_to_raw,geom);
            catch
                disp('Unable to sort!');
            end
        else
            disp('No ephys!');
        end
    end
end
end

%% function sac_analysis
function sac_analysis(path_data_monkey_sorted, sess_list, params, funcs)

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);
        % try to find fhd file
        path_to_fhd = [sess_path, current_rec, filesep, 'raw_data', filesep];
        fhd_file = dir([path_to_fhd, '*.fhd']);
        if isempty(fhd_file)
            fprintf(' -> NO FHD File!\n');
            continue;
        elseif length(fhd_file) > 1
            fprintf(' -> CONFLICT! Two FHD Files for one recording!');
            continue;
        else
            fprintf("\n           Converting");
        end

        fhd_name = fhd_file(1).name;
        fhd_name = fhd_name(1:end-4);

        % Convert if not converted
        if isempty(dir([path_to_fhd, fhd_name, '.mat']))
            fprintf('\n');
            fhd_to_mat(path_to_fhd);
        else
            fprintf(' -> already done!\n');
        end

        try
        % Analyze with Sac_Sorter
            flag_figure = true;
            [SACS_ALL_DATA, TRIALS_DATA, EXPERIMENT_PARAMS] = ...
                funcs.monkey_behavior_sac([path_to_fhd, fhd_name '.mat'], flag_figure, params, funcs);
        catch
            disp('Error!');
            continue;
        end

        path_to_analyzed_eye = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'eye', filesep];

        mkdir(path_to_analyzed_eye);

        % Save _ANALYZED.mat Data to disk
        save([path_to_analyzed_eye fhd_name '_ANALYZED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-v7.3');

        % Save _REDUCED.mat Data to disk
        rmfields_list = {'eye_l_vm_filt', 'eye_l_vy_filt', 'eye_l_vx_filt', 'eye_l_py_filt', 'eye_l_px_filt', ...
            'eye_r_vm_filt', 'eye_r_vy_filt', 'eye_r_vx_filt', 'eye_r_py_filt', 'eye_r_px_filt', ...
            'time', 'time_1K', 'target_visible', 'reward', 'tgt_py', 'tgt_px', 'time_tgt', ...
            'eye_l_vm', 'eye_r_vm', 'eye_l_vy', 'eye_l_vx', 'eye_r_vy', 'eye_r_vx', ...
            'eye_l_py', 'eye_l_px', 'eye_r_py', 'eye_r_px', 'time_eyelink', 'inds_invalid', 'inds_trial'};
        TRIALS_DATA = rmfield(TRIALS_DATA,rmfields_list);
        save([path_to_analyzed_eye fhd_name '_REDUCED.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', '-v7.3');

        if flag_figure
            path_to_analyzed_eye_fig = [sess_path, current_rec, filesep, 'analyzed_figs', filesep, ...
                'behavior_data', filesep, 'eye', filesep];
            % Save Fig
            hFig_ = gcf;
            mkdir(path_to_analyzed_eye_fig)
            saveas(hFig_,[path_to_analyzed_eye_fig fhd_name '_sac_sorter'], 'pdf');
            saveas(hFig_,[path_to_analyzed_eye_fig fhd_name '_sac_sorter'], 'png');
            close(hFig_);
        end
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function Recalibrate
function recalibrate(path_data_monkey_sorted, sess_list, params, funcs)

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);
    flag_used = ones(length(rec_list),1);

    num_used_rec = 1;

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);

        % Load Sac_Sorter results
        path_to_analyzed_eye = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'eye', filesep];
        % load _ANALYZED.mat Data
        file_name = dir([path_to_analyzed_eye '*_ANALYZED.mat']);
        try
            load([path_to_analyzed_eye file_name(1).name], ...
                'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA');
        catch
            flag_used(counter_rec) = 0;
            fprintf('\n No Behavior!\n')
            continue
        end
        flag_used(counter_rec) = 1;
        TRIALS_DATA_M(num_used_rec) = TRIALS_DATA;
        SACS_ALL_DATA_M(num_used_rec) = SACS_ALL_DATA;
        EXPERIMENT_PARAMS_M(num_used_rec) = EXPERIMENT_PARAMS;
        num_used_rec = num_used_rec + 1;
        fprintf('\n load Complete!\n')
    end
    flag_use_avg = 1;
    min_num_sample = 1;

    CAL_MATRIX =  eye(3);
    while(1)
        [SACS_ALL_DATA_M,TRIALS_DATA_M,EXPERIMENT_PARAMS_M, CAL_MATRIX_M] = ...
            funcs.Sac_Recal(SACS_ALL_DATA_M,...
            TRIALS_DATA_M,EXPERIMENT_PARAMS_M,...
            flag_use_avg,min_num_sample, params, funcs);
        CAL_MATRIX = CAL_MATRIX*CAL_MATRIX_M;
        run_recal_resp = questdlg('Re-run calibration?');
        if strcmp(run_recal_resp,'No')
            break;
        elseif strcmp(run_recal_resp,'Cancel')
            return;
        end

    end

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        if flag_used(counter_rec) == 0
            continue;
        end
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);

        % Save Sac_Sorter results
        path_to_analyzed_eye = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'eye', filesep];
        % load _ANALYZED.mat Data
        TRIALS_DATA = TRIALS_DATA_M(sum(flag_used(1:counter_rec)));
        SACS_ALL_DATA = SACS_ALL_DATA_M(sum(flag_used(1:counter_rec)));
        EXPERIMENT_PARAMS = EXPERIMENT_PARAMS_M(sum(flag_used(1:counter_rec)));
        file_name = dir([path_to_analyzed_eye '*_ANALYZED.mat']);
        name_ = file_name(1).name(1:end-4);
        save([path_to_analyzed_eye name_ '_RECAL.mat'], ...
            'EXPERIMENT_PARAMS', 'TRIALS_DATA', 'SACS_ALL_DATA', 'CAL_MATRIX','-v7.3');
        fprintf('\n Save Complete!\n')
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function lick_analysis
function lick_analysis(path_data_monkey_sorted, sess_list, params, funcs)

if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);
       
        % try to find DLC csv file
        path_to_csv = [sess_path, current_rec, filesep, 'analyzed_data',...
            filesep, 'behavior_data',filesep, 'tongue', filesep];
        csv_file = dir([path_to_csv, '*_DLC.csv']);
        if isempty(csv_file)
            fprintf(' -> NO DLC csv File!\n');
            continue;
        elseif length(csv_file) > 1
            fprintf(' -> CONFLICT! Two DLC csv Files for one recording!');
            continue;
        else
            fprintf("\n           Analyzing");
        end

        csv_name = csv_file(1).name;
        [~,csv_name,~] = fileparts(csv_name);   

        % Convert to mat if not converted
        if isempty(dir([path_to_csv, csv_name, '.mat']))
            fprintf('\n');
            PGH_csv_to_mat(path_to_csv);
        else
            fprintf(' -> already done!\n');
        end

        % Analyze with Lick_Sorter
        flag_figure = true;
        [LICKS_ALL_DATA, EXPERIMENT_PARAMS] = ...
            funcs.monkey_behavior_lick([path_to_csv, csv_name '.mat'], flag_figure, params, funcs);

        path_to_analyzed_tongue = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'tongue', filesep];

        mkdir(path_to_analyzed_tongue);

        % Save _ANALYZED.mat Data to disk
        save([path_to_analyzed_tongue csv_name(1:13) '_ANALYZED.mat'], ...
            'EXPERIMENT_PARAMS', 'LICKS_ALL_DATA', '-v7.3');

      

        if flag_figure
            path_to_analyzed_tongue = [sess_path, current_rec, filesep, 'analyzed_figs', filesep, ...
                'behavior_data', filesep, 'tongue', filesep];
            % Save Fig
            hFig_ = gcf;
            mkdir(path_to_analyzed_tongue)
            saveas(hFig_,[path_to_analyzed_tongue csv_name(1:13) '_lick_sorter'], 'pdf');
            saveas(hFig_,[path_to_analyzed_tongue csv_name(1:13) '_lick_sorter'], 'png');
            close(hFig_);
        end
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function make metadata
function generate_session_meta_data(path_data_monkey_sorted, sess_list, params, funcs)

if ~strcmp(path_data_monkey_sorted(end), filesep)
    path_data_monkey_sorted = [path_data_monkey_sorted filesep];
end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    task_num = input('Task (rand_corr: 1, fixed_corr: 2, back_step:3, Xaxis_single: 4): ');
    alignment_ = input('Alignment: ', 's');
    elec_num = input('Electrode (M1: 1, M2: 2, Tetrode: 3, Heptode: 4): ');

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);
    flag_eye = zeros(length(rec_list),1);
    flag_ephys = zeros(length(rec_list),1);
    flag_tongue = zeros(length(rec_list),1);
    trial_number = zeros(length(rec_list),1);
    task = zeros(length(rec_list),1);
    alignment = cell(length(rec_list),1);
    mpm = nan(length(rec_list),1);
    elec = nan(length(rec_list),1);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' rec ', current_rec ' \n']);

        % Check Sac Sorter Results
        path_to_analyzed_eye = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'eye', filesep];
        % check _ANALYZED_RECAL.mat Data
        eye_file = dir([path_to_analyzed_eye '*_ANALYZED_RECAL.mat']);
        flag_eye(counter_rec) = not(isempty(eye_file));
        if flag_eye(counter_rec)
            load([path_to_analyzed_eye eye_file(1).name],'EXPERIMENT_PARAMS');
            % trial number
            trial_number(counter_rec) = EXPERIMENT_PARAMS.num_trials;
            task(counter_rec) = task_num;
        else
            trial_number(counter_rec) = 0;
            task(counter_rec) = nan;
        end

        % Check Units Directory
        path_to_sorted_data = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'sorted_data', filesep];
        % check the units folder content
        flag_ephys(counter_rec) = logical(exist(path_to_sorted_data,'dir'));

        % Check Lick Sorter Results
        path_to_analyzed_lick = [sess_path, current_rec, filesep, 'analyzed_data', filesep, ...
            'behavior_data', filesep, 'tongue', filesep];
        % check _ANALYZED.mat Data
        tongue_file = dir([path_to_analyzed_lick '*_ANALYZED.mat']);
        flag_tongue(counter_rec) = not(isempty(tongue_file));


        if flag_ephys(counter_rec)
            alignment{counter_rec} = alignment_;
            mpm(counter_rec) = input('MPM: ');
            elec(counter_rec) = elec_num;
        end

    end

    meta_data = table(rec_list',flag_ephys,flag_eye,flag_tongue,...
        trial_number,task,alignment,mpm,elec);
    meta_data.Properties.VariableNames = ...
        {'folder_name','ephys','eye','tongue','num_trial','exp',...
        'alignment','MPM','elec'};
    sess_name = regexprep(current_sess(3:end),'-','');
    writetable(meta_data,[sess_path sess_name],'FileType','spreadsheet')

end
fprintf('### ALL DONE. ###\n')
end

%% function unit_extraction
function unit_extraction(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep)
    path_data_monkey_sorted = [path_data_monkey_sorted filesep];
end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ',...
        num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep,...
        current_sess, filesep];

    sess_name = regexprep(current_sess(3:end),'-','');
    sess_meta_data = readtable([sess_path sess_name '.xls']);

    elec_id_ = sess_meta_data.elec(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));
    elec_id = elec_id_(1);

    switch(elec_id)
        case(1)
            MAF_KS_unit_extraction(sess_path);
        case(2)
            MAF_KS_unit_extraction(sess_path);
        case(3)
            MAF_psort_unit_extraction(sess_path);
        case(4)
            MAF_psort_unit_extraction(sess_path);
    end

end
fprintf('### ALL DONE. ###\n')
end

%% Function Sac_alignment
function sac_alignment(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    sess_name = regexprep(current_sess(3:end),'-','');

    sess_meta_data = readtable([sess_path sess_name '.xls']);
    task = sess_meta_data.exp;
    task(isnan(task)) = [];
    task = task(1);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec '\n']);
        path_to_raw = [sess_path, current_rec, filesep, 'raw_data', filesep];
        path_to_save = [sess_path current_rec filesep 'analyzed_data'...
            filesep 'behavior_data' filesep 'eye' filesep];

        if str2double(sess_name) < 191200
            % back_step
            if task == 3
                ESN_event_alignment_v4_rand_corr_auto_pump(...
                    path_to_raw, path_to_save);
            else
                ESN_event_alignment_v4_bypassADC(...
                    path_to_raw, path_to_save);
            end
        else
            ESN_event_alignment_v4_rand_corr_post201911(...
                path_to_raw, path_to_save);
        end
    end
end
fprintf('### ALL DONE. ###\n')
end

%% Function Lick_alignment
function lick_alignment(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    sess_name = regexprep(current_sess(3:end),'-','');

    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.tongue));

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec '\n']);
        path_to_raw = [sess_path, current_rec, filesep, 'raw_data', filesep];
        path_to_save = [sess_path current_rec filesep 'analyzed_data'...
            filesep 'behavior_data' filesep 'tongue' filesep];
        
        rec_name = regexprep(current_rec(3:end),'-','');

        PGH_event_alignment_LED(path_to_raw, path_to_save, rec_name);
    end
end
fprintf('### ALL DONE. ###\n')
end

%% Function sac_neuro_modulation
function sac_neuro_modulation(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    sess_name = regexprep(current_sess(3:end),'-','');

    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);
        path_to_rec = [sess_path, current_rec, filesep];
        path_to_units = [path_to_rec, 'analyzed_data', filesep,...
            'units', filesep];
        unit_list = dir([path_to_units, ...
            regexprep(current_rec(3:end),'-','') '*']);
        unit_list([unit_list.isdir] == 0) = [];
        unit_list = {unit_list.name};
        num_unit = length(unit_list);

        if num_unit == 0
            fprintf('\n No Units!\n')
            continue;
        end

        for counter_unit = 1 : 1: num_unit
            current_unit = unit_list{counter_unit};
            fprintf(['\n          ', num2str(counter_unit), ' / ' num2str(num_unit) ' Analyzing rec ', current_unit '\n']);
            [SACS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = ...
                funcs.add_ephys_sac_sorter(path_to_rec,current_unit, params, funcs);
            % Save _sac Data to disk
            save([path_to_units current_unit filesep current_unit '_sac.mat'], ...
                'SACS_ALL_DATA', 'Neural_Properties', ...
                'EXPERIMENT_PARAMS', '-v7.3');
        end
    end
end
fprintf('### ALL DONE. ###\n')
end

%% Function lick_neuro_modulation
function lick_neuro_modulation(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end
num_sess = length(sess_list);
% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    rec_list = dir([sess_path , current_sess '*']);
    rec_list([rec_list.isdir] == 0) = [];
    rec_list = {rec_list.name};

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);
        path_to_rec = [sess_path, current_rec, filesep];
        path_to_units = [path_to_rec, 'analyzed_data', filesep,...
            'units', filesep];
        unit_list = dir([path_to_units, ...
            regexprep(current_rec(3:end),'-','') '*']);
        unit_list([unit_list.isdir] == 0) = [];
        unit_list = {unit_list.name};
        num_unit = length(unit_list);

        if num_unit == 0
            fprintf('\n No Units!\n')
            continue;
        end

        for counter_unit = 1 : 1: num_unit
            current_unit = unit_list{counter_unit};
            fprintf(['\n          ', num2str(counter_unit), ' / ' num2str(num_unit) ' Analyzing rec ', current_unit '\n']);
            [LICKS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS] = ...
                funcs.add_ephys_lick_sorter(path_to_rec,current_unit, params, funcs);
            % Save _sac Data to disk
            save([path_to_units current_unit filesep current_unit '_lick.mat'], ...
                'LICKS_ALL_DATA', 'Neural_Properties', ...
                'EXPERIMENT_PARAMS', '-v7.3');
        end
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function rec units summary
function generate_rec_units_summary(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];
    sess_name = regexprep(current_sess(3:end),'-','');

    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    has_tongue = sess_meta_data.tongue(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    num_rec = length(rec_list);

    for counter_rec = 1 : 1 : num_rec
        current_rec = rec_list{counter_rec};
        fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' ']);
        path_to_rec = [sess_path, current_rec, filesep];
        path_to_units = [path_to_rec, 'analyzed_data', filesep,...
            'units', filesep];
        rec_name = regexprep(current_rec(3:end),'-','');
        unit_list = dir([path_to_units, ...
            rec_name, '*']);
        unit_list([unit_list.isdir] == 0) = [];
        unit_list = {unit_list.name};
        num_unit = length(unit_list);

        if num_unit == 0
            fprintf('\n No Units!\n')
            continue;
        end

        clearvars units_summary;

        for counter_unit = 1 : 1: num_unit
            current_unit = unit_list{counter_unit};
            fprintf(['\n          ', num2str(counter_unit), ' / ' num2str(num_unit) ' Analyzing rec ', current_unit '\n']);
            load([path_to_units current_unit filesep current_unit '_sac.mat'], ...
                'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');
            if has_tongue(counter_rec)
                load([path_to_units current_unit filesep current_unit '_lick.mat'], ...
                    'LICKS_ALL_DATA');
            else
                LICKS_ALL_DATA = [];
            end
            unit_summary_ = funcs.unit_summary_func(...
                SACS_ALL_DATA, LICKS_ALL_DATA, Neural_Properties, EXPERIMENT_PARAMS, params, funcs);
            unit_summary_.name = current_unit;
            units_summary(counter_unit) = unit_summary_;
        end
        save([path_to_rec rec_name '_units_summary'], 'units_summary', '-v7.3')
    end
end
fprintf('### ALL DONE. ###\n')
end

%% function combine_bundle_files()
function combine_bundle_files(path_data_monkey_sorted, sess_list, params, funcs)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];

    if exist([sess_path 'units'],'dir')
        rmdir([sess_path 'units'],'s');
    end

    sess_name = regexprep(current_sess(3:end),'-','');
    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    num_recs = length(rec_list);

    bndl = cell(num_recs,1);
    type = cell(num_recs,1);
    unit_names = cell(num_recs,1);
    for counter_recs = 1 : 1 : num_recs
        current_rec = rec_list{counter_recs};
        path_to_rec = [sess_path, current_rec, filesep];
        path_to_rec_units = [path_to_rec, 'analyzed_data', filesep,...
            'units', filesep];
        rec_name = regexprep(current_rec(3:end),'-','');
        unit_list = dir([path_to_rec_units, ...
            rec_name, '*']);
        unit_list([unit_list.isdir] == 0) = [];
        unit_list = {unit_list.name};
        num_unit = length(unit_list);

        % read rec_meta_data
        if isfile([path_to_rec rec_name '.xls'])
            rec_meta_data = ...
                readtable([path_to_rec rec_name '.xls']);
            bndl_ = rec_meta_data.cell;
            if  not(num_unit == length(bndl_))
                error('Inconsistent Recording meta-data!')
            end
            idx_valid = not(bndl_ == -1);
            if not(length(bndl_(idx_valid)) == length(unique(bndl_(idx_valid))))
                error(['duplicate cells in rec ', current_rec]);
            end
            if iscell(rec_meta_data.type)
                type_ = rec_meta_data.type;
            else
                type_ = repmat({''},numel(rec_meta_data.type),1);
            end
            name_ = rec_meta_data.unit_name;
            bndl{counter_recs} = bndl_(idx_valid);
            type{counter_recs} = type_(idx_valid);
            unit_names{counter_recs} = name_(idx_valid);

        else
            error('No rec metadata!');
        end
    end

    cell_ids = unique(cell2mat(bndl));
    num_cells = numel(cell_ids);

    for counter_cell = 1 : 1 : num_cells
        current_cell = cell_ids(counter_cell);
        idx_rec_cell = cellfun(@(x) find(x == current_cell),bndl,'UniformOutput',false);
        idx_rec = cell2mat(idx_rec_cell);
        rec_flag = not(cellfun('isempty',idx_rec_cell));
        cell_recs = rec_list(rec_flag);
        type_ = type(rec_flag);
        unit_names_ = unit_names(rec_flag);
        num_cell_recs = numel(cell_recs);
        clearvars data_recordings
        for counter_recs = 1 : 1 : num_cell_recs
            current_rec = cell_recs{counter_recs};
            current_idx_rec = idx_rec(counter_recs);
            path_to_rec = [sess_path, current_rec, filesep];
            path_to_units = [path_to_rec, 'analyzed_data', filesep,...
                'units', filesep];
            rec_name = regexprep(current_rec(3:end),'-','');

            if counter_recs == 1
                cell_name = [rec_name, ...
                    '_' unit_names_{1}{idx_rec(1)}(1:3) ...
                    num2str(current_cell,'%.3i')];
                fprintf(['      ', num2str(counter_cell), ' / '...
                    num2str(num_cells) ' Analyzing cell ',cell_name '\n']);
            end

            cell_type = type_{counter_recs}{current_idx_rec};
            unit_name = unit_names_{counter_recs}{current_idx_rec};

            % saccade files
            path_to_unit_sac = [path_to_units, ...
                rec_name, '_', unit_name, filesep, ...
                rec_name, '_', unit_name, '_sac.mat'];

            % lick files
            path_to_unit_lick = [path_to_units, ...
                rec_name, '_', unit_name, filesep, ...
                rec_name, '_', unit_name, '_lick.mat'];
            
            data_recording = load(path_to_unit_sac, ...
                'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');

            if isfile(path_to_unit_lick)
                data_recording_lick = load(path_to_unit_lick, ...
                    'LICKS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS');
                data_recording.LICKS_ALL_DATA = data_recording_lick.LICKS_ALL_DATA;
                data_recording.EXPERIMENT_PARAMS.lick_tag_list = data_recording_lick.EXPERIMENT_PARAMS.lick_tag_list;
            end

            data_recording.id             = [rec_name, '_', unit_name];
            data_recording.type           = cell_type;
            data_recordings(counter_recs) = data_recording;
        end
        rec_info.rec_flag = rec_flag;
        rec_info.rec_list = rec_list;
        MAF_combine_recs(data_recordings, sess_path, cell_name, rec_info, params, funcs);
    end

end
fprintf('### ALL DONE. ###\n')
end

%% function extract_cell_interaction
function extract_cell_interactions(path_data_monkey_sorted, sess_list)
if ~strcmp(path_data_monkey_sorted(end), filesep);path_data_monkey_sorted = [path_data_monkey_sorted filesep];end

GLOBAL_XPROB_BINSIZE    = 1e-3;
GLOBAL_XPROB_BEFORE     = 5e-2;
GLOBAL_XPROB_AFTER      = 5e-2;

cross_len = (GLOBAL_XPROB_BEFORE + GLOBAL_XPROB_AFTER)/...
    GLOBAL_XPROB_BINSIZE;

num_sess = length(sess_list);

% Loop over sessions
for counter_sess = 1 : 1 : num_sess
    current_sess = sess_list{counter_sess};
    fprintf(['### ' 'Analyzing sess ', current_sess, ' no. ', num2str(counter_sess), ' / ' num2str(num_sess) ' ### \n']);

    sess_path = [path_data_monkey_sorted, current_sess(1:7), filesep, current_sess, filesep];

    sess_name = regexprep(current_sess(3:end),'-','');
    sess_meta_data = readtable([sess_path sess_name '.xls']);

    rec_list = sess_meta_data.folder_name(logical...
        (sess_meta_data.ephys .* sess_meta_data.eye));

    cell_list = dir([sess_path 'units' filesep sess_name '*']);
    cell_list([cell_list.isdir] == 0) = [];
    cell_list = {cell_list.name};

    cell_ids = cell2mat(cellfun(@(x) str2double(x(18:20)),...
        cell_list,'UniformOutput',false));
    num_cells = numel(cell_ids);
    num_recs = length(rec_list);
    simul_cell_recs = zeros(num_cells,num_cells,num_recs,'logical');

    bndl = cell(num_recs,1);
    unit_names = cell(num_recs,1);
    type = cell(num_recs,1);
    for counter_recs = 1 : 1 : num_recs
        current_rec = rec_list{counter_recs};
        path_to_rec = [sess_path, current_rec, filesep];
        path_to_rec_units = [path_to_rec, 'analyzed_data', filesep,...
            'units', filesep];
        rec_name = regexprep(current_rec(3:end),'-','');
        unit_list = dir([path_to_rec_units, ...
            rec_name, '*']);
        unit_list([unit_list.isdir] == 0) = [];
        unit_list = {unit_list.name};
        num_unit = length(unit_list);

        % read rec_meta_data
        if not(isfile([path_to_rec rec_name '.xls']))
            error('No rec metadata!');
        end

        rec_meta_data = ...
            readtable([path_to_rec rec_name '.xls']);
        bndl_ = rec_meta_data.cell;
        if  not(num_unit == length(bndl_))
            error('Inconsistent Recording meta-data!')
        end
        idx_valid = not(bndl_ == -1);
        bndl_ = bndl_(idx_valid);
        if not(length(bndl_) == length(unique(bndl_)))
            error(['duplicate cells in rec ', current_rec]);
        end
        name_ = rec_meta_data.unit_name;
        name_ = name_(idx_valid);
        if iscell(rec_meta_data.type)
            type_ = rec_meta_data.type;
        else
            type_ = repmat({''},numel(rec_meta_data.type),1);
        end
        type_ = type_(idx_valid);
        bndl{counter_recs} = bndl_;
        unit_names{counter_recs} = name_;
        type{counter_recs} = type_;
        if numel(bndl_) > 1
            pairs = nchoosek(bndl_,2);
        else
            pairs = [];
        end
        for counter_pairs = 1:size(pairs,1)
            current_pair = pairs(counter_pairs,:);
            simul_cell_recs(cell_ids == current_pair(1), ...
                cell_ids == current_pair(2),counter_recs) = 1;
            simul_cell_recs(cell_ids == current_pair(2), ...
                cell_ids == current_pair(1),counter_recs) = 1;
        end
    end

    cross_prob = cell(num_cells);
    for counter_cells1 = 1 : 1 : num_cells
        cell1_id = cell_ids(counter_cells1);
        fprintf(['      ', num2str(counter_cells1), ' / ' num2str(num_cells)...
            ' Analyzing cell ', num2str(cell1_id,'%.3i') ' \n']);
        for counter_cells2 = 1 : 1 : num_cells
            cell2_id = cell_ids(counter_cells2);

            shared_recs = simul_cell_recs(counter_cells1,counter_cells2,:);
            % no shared recs
            if sum(shared_recs) == 0
                cross_prob{counter_cells1,counter_cells2} = [];
                continue;
            end

            % 1: ss1xss2, 2: ss1xcs2, 3: cs1xss2, 4: cs1xcs2 
            cross_prob{counter_cells1,counter_cells2} = nan(5,cross_len+1);

            shared_recs_ind = find(shared_recs);
            ss_index1 = [];
            cs_index1 = [];
            ss_index2 = [];
            cs_index2 = [];
            type_cell1 = type{shared_recs_ind(1)}...
                {bndl{shared_recs_ind(1)} == cell1_id};
            type_cell2 = type{shared_recs_ind(1)}...
                {bndl{shared_recs_ind(1)} == cell2_id};
            if strcmp(type_cell1,'PC')
                flag_ss1 = 1;
                flag_cs1 = 1;
            elseif strcmp(type_cell1,'CS')
                flag_cs1 = 1;
                flag_ss1 = 0;
            else
                flag_ss1 = 1;
                flag_cs1 = 0;
            end
            if strcmp(type_cell1,'PC')
                flag_ss2 = 1;
                flag_cs2 = 1;
            elseif strcmp(type_cell2,'CS')
                flag_cs2 = 1;
                flag_ss2 = 0;
            else
                flag_ss2 = 1;
                flag_cs2 = 0;
            end
            % loop over shared recs
            for counter_recs = 1 : 1 : numel(shared_recs_ind)
                current_rec_ind = shared_recs_ind(counter_recs);
                current_rec = rec_list{current_rec_ind};
                rec_name = regexprep(current_rec(3:end),'-','');
                path_to_rec = [sess_path, current_rec, filesep];
                path_to_rec_units = [path_to_rec, 'analyzed_data', filesep,...
                    'units', filesep];
                % loading ss and cs indices
                unit_name1 = unit_names{current_rec_ind}...
                    {bndl{current_rec_ind} == cell1_id};
                unit_name_cell1 = [rec_name '_' unit_name1];
                unit_name2 = unit_names{current_rec_ind}...
                    {bndl{current_rec_ind} == cell2_id};
                unit_name_cell2 = [rec_name '_' unit_name2];

                file_name = dir([path_to_rec_units unit_name_cell1 ...
                    filesep unit_name_cell1 '_sorted_*.mat']);
                load([path_to_rec_units unit_name_cell1 ...
                    filesep file_name(1).name],...
                    'ss_index','cs_index');

                ss_index1 = [ss_index1; ss_index];
                cs_index1 = [cs_index1; cs_index];

                file_name = dir([path_to_rec_units unit_name_cell2 ...
                    filesep unit_name_cell2 '_sorted_*.mat']);
                load([path_to_rec_units unit_name_cell2 ...
                    filesep file_name(1).name],...
                    'ss_index','cs_index','sample_rate');
                ss_index2 = [ss_index2; ss_index];
                cs_index2 = [cs_index2; cs_index];
            end

            if flag_ss1
                if flag_ss2
                    [xprob, xprob_span] = ESN_extract_xprob(ss_index1,...
                        ss_index2, sample_rate, ...
                        GLOBAL_XPROB_BINSIZE, ...
                        GLOBAL_XPROB_BEFORE, ...
                        GLOBAL_XPROB_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    cross_prob{counter_cells1,counter_cells2}(1,:) = ...
                        xprob_mean;
                    cross_prob{counter_cells1,counter_cells2}(5,:) = ...
                        xprob_span_mean;
                end
                if flag_cs2
                    [xprob, xprob_span] = ESN_extract_xprob(ss_index1,...
                        cs_index2, sample_rate, ...
                        GLOBAL_XPROB_BINSIZE, ...
                        GLOBAL_XPROB_BEFORE, ...
                        GLOBAL_XPROB_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    cross_prob{counter_cells1,counter_cells2}(2,:) = ...
                        xprob_mean;
                    cross_prob{counter_cells1,counter_cells2}(5,:) = ...
                        xprob_span_mean;
                end
            end
            if flag_cs1
                if flag_ss2
                    [xprob, xprob_span] = ESN_extract_xprob(cs_index1,...
                        ss_index2, sample_rate, ...
                        GLOBAL_XPROB_BINSIZE, ...
                        GLOBAL_XPROB_BEFORE, ...
                        GLOBAL_XPROB_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    cross_prob{counter_cells1,counter_cells2}(3,:) = ...
                        xprob_mean;
                    cross_prob{counter_cells1,counter_cells2}(5,:) = ...
                        xprob_span_mean;
                end
                if flag_cs2
                    [xprob, xprob_span] = ESN_extract_xprob(cs_index1,...
                        cs_index2, sample_rate, ...
                        GLOBAL_XPROB_BINSIZE, ...
                        GLOBAL_XPROB_BEFORE, ...
                        GLOBAL_XPROB_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    cross_prob{counter_cells1,counter_cells2}(4,:) = ...
                        xprob_mean;
                    cross_prob{counter_cells1,counter_cells2}(5,:) = ...
                        xprob_span_mean;
                end
            end
        end
    end
    save([sess_path 'units' filesep sess_name '.mat'],'cross_prob',...
        'simul_cell_recs','cell_ids')
end
fprintf('### ALL DONE. ###\n')
end
