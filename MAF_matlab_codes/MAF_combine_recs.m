%% function combine_recs
function MAF_combine_recs(data_recordings, sess_path, unit_name, rec_info, params, funcs)

% Init data_cell
data_cell = struct;

num_recording = numel(data_recordings);

% id
data_cell.id = cell(num_recording, 1);
for counter_recording = 1 : 1 : num_recording
    data_cell.id{counter_recording, 1} = data_recordings(counter_recording).id;
end

% EXPERIMENT_PARAMS
field_names_EXPERIMENT_PARAMS = fieldnames(data_recordings(1).EXPERIMENT_PARAMS);
for counter_field = 1 : 1 : length(field_names_EXPERIMENT_PARAMS)
    field_name_EXPERIMENT_PARAMS = field_names_EXPERIMENT_PARAMS{counter_field};
    data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS) = cell(num_recording, 1);
    for counter_recording = 1 : 1 : num_recording
        data_cell.EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS){counter_recording, 1} = data_recordings(counter_recording).EXPERIMENT_PARAMS.(field_name_EXPERIMENT_PARAMS);
    end
end

data_cell.EXPERIMENT_PARAMS.name = unit_name;
data_cell.EXPERIMENT_PARAMS.rec_info = rec_info;

rec_numbers = find(rec_info.rec_flag == 1);

% Neural_Properties
% Init variables
data_cell.Neural_Properties.SS_num = 0;
data_cell.Neural_Properties.SS_duration = 0;
data_cell.Neural_Properties.SS_firing_rate = 0;
data_cell.Neural_Properties.CS_num = 0;
data_cell.Neural_Properties.CS_firing_rate = 0;
data_cell.Neural_Properties.SS_time = [];
data_cell.Neural_Properties.CS_time = [];
variable_names = {'Corr_data_CS_inds_span', 'Corr_data_CS_bin_size_time',...
    'Corr_data_SS_inds_span', 'Corr_data_SS_bin_size_time', ...
    'Corr_data_SS_SSxSS_AUTO','Corr_data_CS_CSxSS_AUTO'};
for counter_variable = 1 : length(variable_names)
    variable_name = variable_names{counter_variable};
    data_cell.Neural_Properties.(variable_name) = ...
        zeros(size(data_recordings(1).Neural_Properties.(variable_name)));
end

data_cell.Neural_Properties.waveform.ch_map = ...
    data_recordings(1).Neural_Properties.waveform.ch_map;

data_cell.Neural_Properties.waveform.ss_wave = ...
    data_recordings(1).Neural_Properties.waveform.ss_wave.*...
    data_recordings(1).Neural_Properties.SS_num;
data_cell.Neural_Properties.waveform.cs_wave = ...
    data_recordings(1).Neural_Properties.waveform.cs_wave.*...
    data_recordings(1).Neural_Properties.CS_num;

data_cell.Neural_Properties.waveform.ss_wave_ch = ...
    data_recordings(1).Neural_Properties.waveform.ss_wave_ch;
data_cell.Neural_Properties.waveform.cs_wave_ch = ...
    data_recordings(1).Neural_Properties.waveform.cs_wave_ch;

% cell type
data_cell.Neural_Properties.type = ...
    data_recordings(1).Neural_Properties.type;

% Loop over recordings
for counter_recording = 1 : 1 : num_recording
    % type
    current_type = data_recordings(counter_recording).type;
    if not(strcmp(current_type,''))
        data_cell.Neural_Properties.type = current_type;
    end
    % waveform
    if counter_recording > 1
        current_ss_wave_ch = data_recordings...
            (counter_recording).Neural_Properties.waveform.ss_wave_ch;
        current_cs_wave_ch = data_recordings...
            (counter_recording).Neural_Properties.waveform.cs_wave_ch;
        current_ss_wave = data_recordings...
            (counter_recording).Neural_Properties.waveform.ss_wave;
        current_cs_wave = data_recordings...
            (counter_recording).Neural_Properties.waveform.cs_wave;

        old_ss_wave_ch = data_cell.Neural_Properties.waveform.ss_wave_ch;
        old_cs_wave_ch = data_cell.Neural_Properties.waveform.cs_wave_ch;
        old_ss_wave = data_cell.Neural_Properties.waveform.ss_wave;
        old_cs_wave = data_cell.Neural_Properties.waveform.cs_wave;

        [new_ss_wave_ch,idx_ss_old,idx_ss_current] = ...
            intersect(old_ss_wave_ch,current_ss_wave_ch);
        [new_cs_wave_ch,idx_cs_old,idx_cs_current] = ...
            intersect(old_cs_wave_ch,current_cs_wave_ch);

        ss_num = data_recordings(...
            counter_recording).Neural_Properties.SS_num;
        cs_num = data_recordings(...
            counter_recording).Neural_Properties.CS_num;

        data_cell.Neural_Properties.waveform.ss_wave_ch = new_ss_wave_ch;
        data_cell.Neural_Properties.waveform.cs_wave_ch = new_cs_wave_ch;
        data_cell.Neural_Properties.waveform.ss_wave = nansum(...
            cat(3,old_ss_wave(idx_ss_old,:), ...
            current_ss_wave(idx_ss_current,:).*ss_num),3);
        data_cell.Neural_Properties.waveform.cs_wave = nansum(...
            cat(3,old_cs_wave(idx_cs_old,:), ...
            current_cs_wave(idx_cs_current,:).*cs_num),3);
    end

    % Corr_data_CS and Corr_data_SS Vars
    data_cell.Neural_Properties.SS_num = data_cell.Neural_Properties.SS_num + ...
        data_recordings(counter_recording).Neural_Properties.SS_num;
    data_cell.Neural_Properties.SS_duration = data_cell.Neural_Properties.SS_duration + ...
        data_recordings(counter_recording).Neural_Properties.SS_duration;
    data_cell.Neural_Properties.CS_num = data_cell.Neural_Properties.CS_num + ...
        data_recordings(counter_recording).Neural_Properties.CS_num;

    data_cell.Neural_Properties.SS_time = vertcat(data_cell.Neural_Properties.SS_time ,...
        data_recordings(counter_recording).Neural_Properties.SS_time);
    data_cell.Neural_Properties.CS_time = vertcat(data_cell.Neural_Properties.CS_time, ...
        data_recordings(counter_recording).Neural_Properties.CS_time);
    % compute weighted average for these variables
    for counter_variable = 1 : length(variable_names)
        variable_name = variable_names{counter_variable};
        if contains(variable_name, 'Corr_data_SS')
            num_spike = data_recordings(counter_recording).Neural_Properties.SS_num;
        elseif contains(variable_name, 'Corr_data_CS')
            num_spike = data_recordings(counter_recording).Neural_Properties.CS_num;
        end
        data_cell.Neural_Properties.(variable_name) = nansum(cat(3,data_cell.Neural_Properties.(variable_name), ...
            ( data_recordings(counter_recording).Neural_Properties.(variable_name) .* num_spike)),3);
    end
end
data_cell.Neural_Properties.SS_firing_rate = ...
    data_cell.Neural_Properties.SS_num ./ data_cell.Neural_Properties.SS_duration;
data_cell.Neural_Properties.CS_firing_rate = ...
    data_cell.Neural_Properties.CS_num ./ data_cell.Neural_Properties.SS_duration;
% devide by number of events to compute the weigted average
for counter_variable = 1 : length(variable_names)
    variable_name = variable_names{counter_variable};
    if contains(variable_name, 'Corr_data_SS')
        num_spike = data_cell.Neural_Properties.SS_num;
    elseif contains(variable_name, 'Corr_data_CS')
        num_spike = data_cell.Neural_Properties.CS_num;
    end
    data_cell.Neural_Properties.(variable_name) = ...
        data_cell.Neural_Properties.(variable_name) ./ num_spike;
end

data_cell.Neural_Properties.waveform.ss_wave = ...
    data_cell.Neural_Properties.waveform.ss_wave./data_cell.Neural_Properties.SS_num;
data_cell.Neural_Properties.waveform.cs_wave = ...
    data_cell.Neural_Properties.waveform.cs_wave./data_cell.Neural_Properties.CS_num;

% SACS_ALL_DATA
field_names_SACS_ALL_DATA = fieldnames(data_recordings(1).SACS_ALL_DATA);
for counter_field = 1 : 1 : length(field_names_SACS_ALL_DATA)
    field_name_SACS_ALL_DATA = field_names_SACS_ALL_DATA{counter_field};
    data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA) = horzcat(...
            data_cell.SACS_ALL_DATA.(field_name_SACS_ALL_DATA), ...
            data_recordings(counter_recording).SACS_ALL_DATA.(field_name_SACS_ALL_DATA));
    end
end

data_cell.SACS_ALL_DATA.rec_num = [];
for counter_recording = 1 : 1 : num_recording
    data_cell.SACS_ALL_DATA.rec_num = horzcat(...
        data_cell.SACS_ALL_DATA.rec_num, ...
        rec_numbers(counter_recording)*...
        ones(size(data_recordings(counter_recording).SACS_ALL_DATA.validity)));
end

% LICKS_ALL_DATA
if isfield(data_recordings(1),'LICKS_ALL_DATA')
    field_names_LICKS_ALL_DATA = fieldnames(data_recordings(1).LICKS_ALL_DATA);
    for counter_field = 1 : 1 : length(field_names_LICKS_ALL_DATA)
        field_name_LICKS_ALL_DATA = field_names_LICKS_ALL_DATA{counter_field};
        data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = [];
        for counter_recording = 1 : 1 : num_recording
            data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA) = horzcat(...
                data_cell.LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA), ...
                data_recordings(counter_recording).LICKS_ALL_DATA.(field_name_LICKS_ALL_DATA));
        end
    end

    data_cell.LICKS_ALL_DATA.rec_num = [];
    for counter_recording = 1 : 1 : num_recording
        data_cell.LICKS_ALL_DATA.rec_num = horzcat(...
            data_cell.LICKS_ALL_DATA.rec_num, ...
            rec_numbers(counter_recording)*...
            ones(size(data_recordings(counter_recording).LICKS_ALL_DATA.validity)));
    end
    LICKS_ALL_DATA = data_cell.LICKS_ALL_DATA;
end

% Save data_cell
id = data_cell.id;
EXPERIMENT_PARAMS = data_cell.EXPERIMENT_PARAMS;
Neural_Properties = data_cell.Neural_Properties;
SACS_ALL_DATA = data_cell.SACS_ALL_DATA;
cell_name = [unit_name '_' 'combine' '_' num2str(length(data_cell.id))];
fprintf([cell_name ': Saving .mat file ...']);
save_path = [sess_path 'units' filesep cell_name filesep];
mkdir(save_path);
if exist('LICKS_ALL_DATA','var')
    save([save_path cell_name '.mat'],...
        'SACS_ALL_DATA', 'LICKS_ALL_DATA', 'Neural_Properties',...
        'EXPERIMENT_PARAMS', 'id','-v7.3')
else
    save([save_path cell_name '.mat'],...
        'SACS_ALL_DATA', 'Neural_Properties', 'EXPERIMENT_PARAMS',...
        'id','-v7.3')
end
fprintf(' --> Completed. \n')

% plot Neural_Properties
params.cell_name    = cell_name;
params.duration     = Neural_Properties.SS_duration;
params.num_trials   = nansum(cell2mat(EXPERIMENT_PARAMS.num_trials));
fprintf([cell_name ': Saving neural properties plot ...'])
MAF_plot_NeuralProperties(Neural_Properties, params)
hFig_ = gcf;
figs_path = [sess_path 'units' filesep cell_name filesep...
    'analyzed_figs' filesep];
mkdir(figs_path);
saveas(hFig_,[figs_path cell_name '_Neural_Properties.pdf'], 'pdf');
close(hFig_)
fprintf(' --> Completed. \n')

% plot_sac_sorter
params.cell_name    = cell_name;
params.duration     = Neural_Properties.SS_duration;
params.sac_tag_list = EXPERIMENT_PARAMS.sac_tag_list{1};
params.num_trials   = nansum(cell2mat(EXPERIMENT_PARAMS.num_trials));
fprintf([cell_name ': Saving saccade sorter plot ...'])
ESN_plot_sac_sorter(SACS_ALL_DATA, params)
hFig_ = gcf;
figs_path = [sess_path 'units' filesep cell_name filesep...
    'analyzed_figs' filesep 'behavior_data' filesep 'eye' filesep];
mkdir(figs_path);
saveas(hFig_,[figs_path cell_name '_sac_sorter'], 'pdf');
close(hFig_)
fprintf(' --> Completed. \n')

% plot_lick_sorter
if exist('LICKS_ALL_DATA','var')
    params.cell_name    = cell_name;
    params.duration     = Neural_Properties.SS_duration;
    params.lick_tag_list = EXPERIMENT_PARAMS.lick_tag_list{1};
    fprintf([cell_name ': Saving lick sorter plot ...'])
    PGH_plot_lick_sorter(LICKS_ALL_DATA, params)
    hFig_ = gcf;
    figs_path = [sess_path 'units' filesep cell_name filesep...
        'analyzed_figs' filesep 'behavior_data' filesep 'tongue' filesep];
    mkdir(figs_path);
    saveas(hFig_,[figs_path cell_name '_lick_sorter'], 'pdf');
    close(hFig_)
    fprintf(' --> Completed. \n')
end

end
