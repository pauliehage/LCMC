%% function ESN_plot_neural_modulation_v4(num_data_set)
function PH_plot_DCN_modulation(num_data_set, params)
if nargin < 1
    num_data_set = 1;
    params.manual = 1;
    params.auto = 1;
elseif nargin == 1
    params.manual = 1;
    params.auto = 1;
elseif nargin > 1
    params.auto = 1;
    params.manual = 0;
    
end
%% Build EPHYS_ and BEHAVE_ for each single dataset
clearvars EPHYS BEHAVE
for counter_dataset = 1 : 1 : num_data_set
    [EPHYS_(counter_dataset), BEHAVE_(counter_dataset)] = build_EPHYS_BEHAVE_single_dataset(num_data_set,params);
end
EPHYS.CH_sorted_file_name = EPHYS_(1).CH_sorted_file_name;
EPHYS.CH_sorted_file_path = EPHYS_(1).CH_sorted_file_path;

[~, file_name, ~]  = fileparts(EPHYS.CH_sorted_file_name);
if num_data_set > 1
    file_name = [file_name '_combine_' num2str(num_data_set)];
end
Neural_Properties_data.file_name = file_name;
Neural_Properties_data.file_path = EPHYS.CH_sorted_file_path;

%% Plot-1 SS train cue_present 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_cue_present_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_cue_present_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'prim');
end

raster_data_cue_present = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_cue_present.fig_num_               = 1;
plot_data_cue_present.xlabel_text_raster_    = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
plot_data_cue_present.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data_cue_present.inds_span = EPHYS_(1).CH_EVE.inds_span_cue_present;
fig_handle_(plot_data_cue_present.fig_num_)  = plot_rasters_data(raster_data_cue_present, plot_data_cue_present);

sgtitle(fig_handle_(plot_data_cue_present.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-2 SS & CS train primSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_primSac_onset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_primSac_onset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'prim');
end

raster_data_primSac_onset  = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_primSac_onset.fig_num_               = 2;
plot_data_primSac_onset.xlabel_text_raster_    = {'Time relative to prim sac onset (ms)', 'Directions based on primary sac'};
plot_data_primSac_onset.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data_primSac_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_onset;
fig_handle_(plot_data_primSac_onset.fig_num_)  = plot_rasters_data(raster_data_primSac_onset , plot_data_primSac_onset);

sgtitle(fig_handle_(plot_data_primSac_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-3 SS & CS train primSac_offset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_primSac_offset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_primSac_offset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'corr');
end

raster_data_primSac_offset = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_primSac_offset.fig_num_               = 3;
plot_data_primSac_offset.xlabel_text_raster_    = {'Time relative to prim sac offset (ms)', 'Directions based on secondary sac'};
plot_data_primSac_offset.xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
plot_data_primSac_offset.inds_span = EPHYS_(1).CH_EVE.inds_span_primSac_offset;
fig_handle_(plot_data_primSac_offset.fig_num_)  = plot_rasters_data(raster_data_primSac_offset, plot_data_primSac_offset);

sgtitle(fig_handle_(plot_data_primSac_offset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-4 SS & CS train corrSac_onset 
clearvars raster_data_
for counter_dataset = 1 : 1 : num_data_set
    EPHYS_inds_event            = EPHYS_(counter_dataset).CH_EVE.EPHYS_EB_aligned_inds_corrSac_onset_1K;
    BEHAVE_inds_event           = EPHYS_(counter_dataset).CH_EVE.BEHAVE_EB_aligned_inds_corrSac_onset_1K;
    raster_data_(counter_dataset).data = single_dataset_raster(...
        EPHYS_(counter_dataset), BEHAVE_(counter_dataset), ...
        EPHYS_inds_event, BEHAVE_inds_event, 'corr');
end

raster_data_corrSac_onset = concatenate_dataset(raster_data_, 'data', @vertcat);

plot_data_corrSac_onset.fig_num_               = 4;
plot_data_corrSac_onset.xlabel_text_raster_    = {'Time relative to corr sac onset (ms)', 'Directions based on secondary sac'};
plot_data_corrSac_onset.xlabel_text_CS_probab_ = {'CS probability based on [-200 0]ms'};
plot_data_corrSac_onset.inds_span = EPHYS_(1).CH_EVE.inds_span_corrSac_onset;
fig_handle_(plot_data_corrSac_onset.fig_num_)  = plot_rasters_data(raster_data_corrSac_onset, plot_data_corrSac_onset);

sgtitle(fig_handle_(plot_data_corrSac_onset.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Plot-5 Neural Properties
CH_sorted_ = concatenate_dataset(EPHYS_, 'CH_sorted', @horzcat);
EPHYS.CH_sorted.SS_data   = concatenate_dataset(CH_sorted_.SS_data, [], @vertcat);
EPHYS.CH_sorted.CS_data   = concatenate_dataset(CH_sorted_.CS_data, [], @vertcat);
EPHYS.CH_sorted.Corr_data = concatenate_dataset(CH_sorted_.Corr_data, [], @vertcat);

if isfield(EPHYS.CH_sorted.SS_data, 'SS_waveform_hipass')
    EPHYS.CH_sorted.SS_data.SS_waveform = EPHYS.CH_sorted.SS_data.SS_waveform_hipass;
    EPHYS.CH_sorted.CS_data.CS_waveform = EPHYS.CH_sorted.CS_data.CS_waveform_hipass;
end

Neural_Properties_data.CH_sorted.SS_data   = EPHYS.CH_sorted.SS_data;
Neural_Properties_data.CH_sorted.CS_data   = EPHYS.CH_sorted.CS_data;
Neural_Properties_data.CH_sorted.Corr_data = EPHYS.CH_sorted.Corr_data;

% figure
Line_Color = lines(7);
plot_data.fig_num_ = 5;
fig_handle_(plot_data.fig_num_) = figure(plot_data.fig_num_);
clf(fig_handle_(plot_data.fig_num_))

% subplot(2,2,1) Waveform
plot_handle_(1) = subplot(2,2,1);
hold on

plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)+std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform)-std(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(1,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.SS_data.SS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(1,:))

plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)+std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform)-std(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 1, 'Color', Line_Color(7,:))
plot((1:180)/30-2, mean(Neural_Properties_data.CH_sorted.CS_data.CS_waveform), ...
    'LineWidth', 2, 'Color', Line_Color(7,:))

xlabel('Time (ms)')
ylabel('Voltage (uv)')
title('CS & SS Waveform')

% subplot(2,2,2) Probablities
plot_handle_(2) = subplot(2,2,2);
hold on

SSxSS_AUTO = Neural_Properties_data.CH_sorted.Corr_data.SS_SSxSS_AUTO;
if size(SSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span);
else
    inds_span =      Neural_Properties_data.CH_sorted.Corr_data.SS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.SS_bin_size_time);
if (~isempty(SSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(SSxSS_AUTO, 1) > 1
        prob_value_ = mean(SSxSS_AUTO);
    else
        prob_value_ =      SSxSS_AUTO;
    end
    prob_value_(round(length(inds_span)/2)) = 0;
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(SSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(SSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(SSxSS_AUTO));
    y_axis_mean_ = nan(size(SSxSS_AUTO));
    y_axis_stdv_ = nan(size(SSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(1,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(1,:))

CSxSS_AUTO    = Neural_Properties_data.CH_sorted.Corr_data.CS_CSxSS_AUTO;
if size(CSxSS_AUTO, 1) > 1
    inds_span = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span);
else
    inds_span      = Neural_Properties_data.CH_sorted.Corr_data.CS_inds_span;
end
bin_size_time = mean(Neural_Properties_data.CH_sorted.Corr_data.CS_bin_size_time);
if (~isempty(CSxSS_AUTO))
    x_axis_ = (inds_span * bin_size_time * 1000)';
    if size(CSxSS_AUTO, 1) > 1
        prob_value_ = mean(CSxSS_AUTO);
    else
        prob_value_ =      CSxSS_AUTO;
    end
    y_axis_mean_ = prob_value_;
    y_axis_stdv_ = sqrt(size(CSxSS_AUTO, 1) .* ( prob_value_ .* (1 - prob_value_) ) ) ./ size(CSxSS_AUTO, 1);
    xlim([x_axis_(1) x_axis_(end)]);
else
    x_axis_      = nan(size(CSxSS_AUTO));
    y_axis_mean_ = nan(size(CSxSS_AUTO));
    y_axis_stdv_ = nan(size(CSxSS_AUTO));
end
plot(x_axis_, y_axis_mean_+y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_-y_axis_stdv_, 'LineWidth', 1, 'Color', Line_Color(7,:))
plot(x_axis_, y_axis_mean_,              'LineWidth', 2, 'Color', Line_Color(7,:))

ylim([0 inf])
xlabel('Time (ms)')
ylabel('Probability')
title('X Probability')

% subplot(2,2,3) ISI SS
plot_handle_(3) = subplot(2,2,3);
hold on

edges_SS = (0 : 0.002 : 0.050) *1000;
ISI_SS = diff(Neural_Properties_data.CH_sorted.SS_data.SS_time) * 1000;
histogram(ISI_SS,edges_SS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(1,:));
set(plot_handle_(3), 'XTick', [0 0.025 0.050]*1000)
xlabel('Time (ms)')
ylabel('Probability')
title('SS Inter-Spike-Interval')

% subplot(2,2,4) ISI CS
plot_handle_(4) = subplot(2,2,4);
hold on
edges_CS = 0 : 0.200 : 5.000;
ISI_CS = diff(Neural_Properties_data.CH_sorted.CS_data.CS_time);
histogram(ISI_CS,edges_CS, 'DisplayStyle', 'stairs', 'Normalization', 'probability', 'LineWidth', 2, 'EdgeColor', Line_Color(7,:));
set(plot_handle_(4), 'XTick', [0 2.5 5.0])
xlabel('Time (s)')
ylabel('Probability')
title('CS Inter-Spike-Interval')

ESN_Beautify_Plot(fig_handle_(plot_data.fig_num_), [8.0 8.0])

sgtitle(fig_handle_(plot_data.fig_num_), Neural_Properties_data.file_name, 'Interpreter', 'none');

%% Save Fig
if params.manual
    response_save_fig = questdlg('Do you want to save the figures?',...
        'Question Dialog','Yes','No','Yes');
else
    response_save_fig = 'No';
    close('all')
end
if contains(response_save_fig, 'Yes')
    fprintf(['Saving plots', ' ...'])
    save_file_path = uigetdir(EPHYS.CH_sorted_file_path, 'Select where to save the figures.');
    if ~isequal(save_file_path,0)
        % saveas(fig_handle_(1),[save_file_path filesep file_name '_modulation_cue_present'], 'pdf');
        saveas(fig_handle_(1),[save_file_path filesep file_name '_modulation_cue_present'], 'png');
        % saveas(fig_handle_(2),[save_file_path filesep file_name '_modulation_primSac_onset'], 'pdf');
        saveas(fig_handle_(2),[save_file_path filesep file_name '_modulation_primSac_onset'], 'png');
        % saveas(fig_handle_(3),[save_file_path filesep file_name '_modulation_primSac_offset'], 'pdf');
        saveas(fig_handle_(3),[save_file_path filesep file_name '_modulation_primSac_offset'], 'png');
        % saveas(fig_handle_(4),[save_file_path filesep file_name '_modulation_corrSac_onset'], 'pdf');
        saveas(fig_handle_(4),[save_file_path filesep file_name '_modulation_corrSac_onset'], 'png');
        saveas(fig_handle_(5),[save_file_path filesep file_name '_properties'], 'pdf');
        saveas(fig_handle_(5),[save_file_path filesep file_name '_properties'], 'png');
        % close(fig_handle_(1))
        % close(fig_handle_(2))
        % close(fig_handle_(3))
        % close(fig_handle_(4))
        % close(fig_handle_(5))
    end
    fprintf(' --> Completed. \n')
end % if contains(response_save_fig, 'Yes')

%% Report Properties
for counter_dataset = 1 : 1 : num_data_set
    [~, file_name, ~]  = fileparts(EPHYS_(counter_dataset).CH_sorted_file_name);
    duration = (EPHYS_(counter_dataset).CH_EVE.EPHYS_time_30K(end)-EPHYS_(counter_dataset).CH_EVE.EPHYS_time_30K(1));
    numCS = length(EPHYS_(counter_dataset).CH_sorted.CS_data.CS_ind);
    freqCS = numCS/duration;
    numSS = length(EPHYS_(counter_dataset).CH_sorted.SS_data.SS_ind);
    freqSS = numSS/duration;
    numTrial = length(BEHAVE_(counter_dataset).TRIALS_DATA.time_end);
    fprintf(['*******************************************' '\n'])
    fprintf([file_name '\n'])
    fprintf([       'dur'   '\t'        'numCS'   '\t'        'freqCS'   '\t'        'numSS'   '\t'        'freqSS'   '\t'        'numTrial'   '\n'])
    fprintf([num2str(duration/60,'%.1f') '\t'  num2str(numCS,'%.0f') '\t' num2str(freqCS,'%.2f') '\t' num2str(numSS,'%.0f') '\t' num2str(freqSS,'%.2f') '\t' num2str(numTrial,'%.0f') '\n'])
end

%% Save plot data
if params.manual
    response_save_data = questdlg('Do you want to save the plot_data?',...
        'Question Dialog','Yes','No','Yes');
else
    response_save_data = 'Yes';
end
if contains(response_save_data, 'Yes')
file_name = [Neural_Properties_data.file_name '_plot_data.mat'];
file_path = Neural_Properties_data.file_path;
if ~params.auto
[save_file_name,save_file_path] = uiputfile([file_path filesep file_name], 'Select where to save the plot data.');
else
save_file_name = file_name;
save_file_path = [file_path filesep '..' filesep 'analyzed_figs'];
end
fprintf(['Saving ' save_file_name ' ... ']);
if ~isequal(save_file_name,0)
save([save_file_path filesep save_file_name], 'Neural_Properties_data', ...
    'raster_data_cue_present',    'plot_data_cue_present',...
    'raster_data_primSac_onset',  'plot_data_primSac_onset', ...
    'raster_data_primSac_offset', 'plot_data_primSac_offset' , ...
    'raster_data_corrSac_onset',  'plot_data_corrSac_onset', ...
    '-v7.3');
end
fprintf(' --> Completed. \n');
end % if contains(response_save_data, 'Yes')

end

%% function ESN_raster_plot_axes
function [x_axis, y_axis] = ESN_raster_plot_axes(train_data_logic, x_axis_values, line_half_len)
if nargin < 2
    x_axis_values = 1 : 1 : size(train_data_logic, 2);
    line_half_len = 0.5;
end
if nargin < 3
    line_half_len = 0.5;
end
train_data_logic = train_data_logic > 0.1;
train_data_row_number = nan(size(train_data_logic));
for counter_row = 1 : size(train_data_logic, 1)
    train_data_row_number(counter_row, train_data_logic(counter_row,:)) = counter_row;
end
train_data_col_number = repmat(x_axis_values(:)', size(train_data_logic,1), 1);
x_axis = [train_data_col_number(:)'; train_data_col_number(:)'; nan(length(train_data_col_number(:)), 1)'];
y_axis = [(train_data_row_number(:)-line_half_len)'; (train_data_row_number(:)+line_half_len)'; nan(length(train_data_row_number(:)), 1)'];
x_axis = x_axis(:);
y_axis = y_axis(:);
end

%% function ESN_smooth
function smooth_data_ = ESN_smooth(data_)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
smooth_data_ = smooth(data_, 21, 'sgolay', 2);
end

%% function build_EPHYS_BEHAVE_single_dataset
function [EPHYS, BEHAVE] = build_EPHYS_BEHAVE_single_dataset(num_data_set,params)

%% load EPHYS sorted DATA
if ~params.auto && num_data_set == 1
    file_path = [pwd filesep];
    [file_name, file_path] = uigetfile([file_path '.psort'], 'Select psort file');
elseif params.auto && num_data_set > 1
    [params.file_name, params.file_path] = uigetfile('.psort', 'Select psort file');
    file_name = params.file_name;
    file_path = params.file_path;
elseif params.auto && ~params.manual && num_data_set == 1
    file_name = [params.file_name '.psort'];
    file_path = params.file_path;
elseif params.auto && params.manual && num_data_set == 1
    [params.file_name, params.file_path] = uigetfile('.psort', 'Select psort file');
    file_name = params.file_name;
    file_path = params.file_path;
    
end
fprintf(['Loading ', file_name, ' ... ']);
% EPHYS.CH_sorted = load([file_path filesep file_name], 'CS_data', 'SS_data');
DATA_PSORT = Psort_read_psort([file_path file_name]);
EPHYS.CH_sorted_file_name = file_name;
EPHYS.CH_sorted_file_path = file_path;
fprintf(' --> Completed. \n')

%% load EPHYS EVENT DATA
if ~params.auto
%     file_name = EPHYS.CH_sorted_file_name(1:13);
    [file_name,file_path] = uigetfile([file_path file_name '_EVE1_aligned.mat'], 'Select EVENT DATA file');
else
    file_name = [EPHYS.CH_sorted_file_name(1:13) '_EVE1_aligned.mat'];
    file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
EPHYS.CH_EVE = load([file_path file_name]);
if isfield(EPHYS.CH_EVE, 'EPHYS_time_15K')
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_15K(:);
else
    EPHYS.CH_EVE.EPHYS_time_30K = EPHYS.CH_EVE.EPHYS_time_30K(:);
end
EPHYS.CH_EVE.EPHYS_time_1K  = EPHYS.CH_EVE.EPHYS_time_1K(:);
EPHYS.CH_EVE.BEHAVE_time_1K = EPHYS.CH_EVE.BEHAVE_time_1K(:);
fprintf(' --> Completed. \n')

%% load BEHAVE DATA
if ~params.auto
% file_name = EPHYS.CH_sorted_file_name(1:13);
[file_name,file_path] = uigetfile([file_path file_name '_ANALYZED.mat'], 'Select _ANALYZED file');
else
file_name = [EPHYS.CH_sorted_file_name(1:13) '_ANALYZED.mat'];
file_path = params.file_path;
end
fprintf(['Loading ', file_name, ' ... ']);
BEHAVE = load([file_path file_name]);
fprintf(' --> Completed. \n')

%% build EPHYS.CH_sorted from DATA_PSORT
ch_data = double(DATA_PSORT.topLevel_data.ch_data);
ch_time = double(DATA_PSORT.topLevel_data.ch_time);
SS_index = find(logical(double(DATA_PSORT.topLevel_data.ss_index)));
CS_index = find(logical(double(DATA_PSORT.topLevel_data.cs_index)));
SS_time = ch_time(SS_index);
CS_time = ch_time(CS_index);

waveform_inds_span = ((-60+1) : 1 : (120));
SS_inds = repmat(waveform_inds_span(:)', length(SS_index), 1) + repmat(SS_index(:), 1, length(waveform_inds_span));
SS_inds(SS_inds < 1) = 1;
SS_inds(SS_inds > length(ch_data)) = length(ch_data);
CS_inds = repmat(waveform_inds_span(:)', length(CS_index), 1) + repmat(CS_index(:), 1, length(waveform_inds_span));
CS_inds(CS_inds < 1) = 1;
CS_inds(CS_inds > length(ch_data)) = length(ch_data);
SS_waveform = ch_data(SS_inds);
CS_waveform = ch_data(CS_inds);

if length(SS_index) == 1
    SS_waveform = SS_waveform(:)';
end

if length(CS_index) == 1
    CS_waveform = CS_waveform(:)';
end

EPHYS.CH_sorted.SS_data.SS_ind = SS_index;
EPHYS.CH_sorted.CS_data.CS_ind = CS_index;
EPHYS.CH_sorted.SS_data.SS_time = SS_time;
EPHYS.CH_sorted.CS_data.CS_time = CS_time;
EPHYS.CH_sorted.SS_data.SS_waveform = SS_waveform;
EPHYS.CH_sorted.CS_data.CS_waveform = CS_waveform;

%% build SSxSS AUTO PROBABILITY
clearvars -except EPHYS BEHAVE
fprintf(['Building SSxSS_AUTO & CSxSS_AUTO PROBABILITY ' ' ...'])
SS_time   = EPHYS.CH_sorted.SS_data.SS_time;
CS_time   = EPHYS.CH_sorted.CS_data.CS_time;
Corr_data = ESN_correlogram(SS_time, CS_time);
EPHYS.CH_sorted.Corr_data.CS_inds_span     = Corr_data.CS_inds_span;
EPHYS.CH_sorted.Corr_data.CS_bin_size_time = Corr_data.CS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_inds_span     = Corr_data.SS_inds_span;
EPHYS.CH_sorted.Corr_data.SS_bin_size_time = Corr_data.SS_bin_size_time;
EPHYS.CH_sorted.Corr_data.SS_SSxSS_AUTO    = Corr_data.SS_SSxSS_AUTO;
EPHYS.CH_sorted.Corr_data.CS_CSxSS_AUTO    = Corr_data.CS_CSxSS_AUTO;

fprintf(' --> Completed. \n')

%% SS & CS train_aligned
clearvars -except EPHYS BEHAVE
fprintf(['Building CS & SS train_aligned', ' ... ']);
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_ = length(EPHYS_time_1K);
CS_time = EPHYS.CH_sorted.CS_data.CS_time;
if isempty(CS_time)
    CS_time = EPHYS_time_1K(1);
end
CS_time(end+1) = max([EPHYS_time_1K(end), CS_time(end)])+1;
SS_time = EPHYS.CH_sorted.SS_data.SS_time;
if isempty(SS_time)
    SS_time = EPHYS_time_1K(1);
end
SS_time(end+1) = max([EPHYS_time_1K(end), SS_time(end)])+1;
EPHYS_CS_train_1K = false(size(EPHYS_time_1K));
EPHYS_SS_train_1K = false(size(EPHYS_time_1K));
counter_CS = find(CS_time >= EPHYS_time_1K(1), 1, 'first');
counter_SS = find(SS_time >= EPHYS_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = EPHYS_time_1K(counter_time_point);
    if time_ponit_>=CS_time(counter_CS)
        EPHYS_CS_train_1K(counter_time_point) = true;
        counter_CS = counter_CS + 1;
    end
    if time_ponit_>=SS_time(counter_SS)
        EPHYS_SS_train_1K(counter_time_point) = true;
        counter_SS = counter_SS + 1;
    end
end
EPHYS.CH_EVE.EPHYS_CS_train_1K = EPHYS_CS_train_1K;
EPHYS.CH_EVE.EPHYS_SS_train_1K = EPHYS_SS_train_1K;

fprintf(' --> Completed. \n')

%% Build BEHAVE_eye_r_vm_filt
clearvars -except EPHYS BEHAVE
fprintf(['Building BEHAVE_eye_r_vm_filt', ' ... ']);
eye_r_vm_filt = cell2mat(BEHAVE.TRIALS_DATA.eye_r_vm_filt(:));
time_1K_cell2mat = cell2mat(BEHAVE.TRIALS_DATA.time_1K(:));
BEHAVE_time_1K     = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_ = length(BEHAVE_time_1K);
BEHAVE_eye_r_vm_filt_1K = nan(size(BEHAVE_time_1K));
time_1K_cell2mat(end+1)    = max([BEHAVE_time_1K(end), time_1K_cell2mat(end)])+1;
counter_time_1K_cell2mat = find(time_1K_cell2mat >= BEHAVE_time_1K(1), 1, 'first');
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_time_1K(counter_time_point);
    if time_ponit_>=time_1K_cell2mat(counter_time_1K_cell2mat)
        BEHAVE_eye_r_vm_filt_1K(counter_time_point) = eye_r_vm_filt(counter_time_1K_cell2mat);
        counter_time_1K_cell2mat = counter_time_1K_cell2mat + 1;
    end
end
EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_1K = BEHAVE_eye_r_vm_filt_1K;
fprintf(' --> Completed. \n')

%% Build Raster Data (inds_span)
clearvars -except EPHYS BEHAVE
fprintf(['Building Raster Plot Data', ' ... ']);
num_trials = length(BEHAVE.TRIALS_DATA.time_end);
BEHAVE_EB_xcorr_time_1K     = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_time_1K;
length_time_ = length(BEHAVE_EB_xcorr_time_1K);
BEHAVE_EB_xcorr_time_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_vmax  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_vmax  = nan(num_trials, 1);
BEHAVE_EB_xcorr_time_corrSac_offset  = nan(num_trials, 1);
for counter_trial = 1 : 1 : num_trials
    BEHAVE_EB_xcorr_time_cue_present(counter_trial) = BEHAVE.TRIALS_DATA.time_state_cue_present{1,counter_trial}(end);
    
    ind_primSac_onset_  = (BEHAVE.SACS_PRIM_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_onset_) ~= 1
        ind_primSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_primSac_onset(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_onset_, counter_trial);
    
    ind_primSac_vmax_  = (BEHAVE.SACS_PRIM_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_vmax_) ~= 1
        ind_primSac_vmax_ = 60;
    end
    BEHAVE_EB_xcorr_time_primSac_vmax(counter_trial)  = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_vmax_, counter_trial);
    
    ind_primSac_offset_ = (BEHAVE.SACS_PRIM_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_PRIM_DATA.inds(:,counter_trial));
    if sum(ind_primSac_offset_) ~= 1
        ind_primSac_offset_ = 150;
    end
    BEHAVE_EB_xcorr_time_primSac_offset(counter_trial) = BEHAVE.SACS_PRIM_DATA.time(ind_primSac_offset_, counter_trial);
    
    ind_corrSac_onset_  = (BEHAVE.SACS_CORR_DATA.ind_start(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_onset_) ~= 1
        ind_corrSac_onset_ = 1;
    end
    BEHAVE_EB_xcorr_time_corrSac_onset(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_onset_, counter_trial);
    
    ind_corrSac_vmax_  = (BEHAVE.SACS_CORR_DATA.ind_vmax(1,counter_trial)  == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_vmax_) ~= 1
        ind_corrSac_vmax_ = 60;
    end
    BEHAVE_EB_xcorr_time_corrSac_vmax(counter_trial)  = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_vmax_, counter_trial);
    
    ind_corrSac_offset_ = (BEHAVE.SACS_CORR_DATA.ind_finish(1,counter_trial) == BEHAVE.SACS_CORR_DATA.inds(:,counter_trial));
    if sum(ind_corrSac_offset_) ~= 1
        ind_corrSac_offset_ = 150;
    end
    BEHAVE_EB_xcorr_time_corrSac_offset(counter_trial) = BEHAVE.SACS_CORR_DATA.time(ind_corrSac_offset_, counter_trial);
end

% num_trials = sum(BEHAVE_EB_xcorr_time_corrSac_offset < BEHAVE_EB_xcorr_time_1K(end));
last_trial_num = find(BEHAVE_EB_xcorr_time_corrSac_offset < BEHAVE_EB_xcorr_time_1K(end), 1, 'last');
first_trial_num = find(BEHAVE_EB_xcorr_time_cue_present > BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
range_trials = first_trial_num : last_trial_num;
num_trials = length(range_trials);

BEHAVE_EB_xcorr_time_cue_present    = BEHAVE_EB_xcorr_time_cue_present(   range_trials);
BEHAVE_EB_xcorr_time_primSac_onset  = BEHAVE_EB_xcorr_time_primSac_onset( range_trials);
BEHAVE_EB_xcorr_time_primSac_vmax   = BEHAVE_EB_xcorr_time_primSac_vmax(  range_trials);
BEHAVE_EB_xcorr_time_primSac_offset = BEHAVE_EB_xcorr_time_primSac_offset(range_trials);
BEHAVE_EB_xcorr_time_corrSac_onset  = BEHAVE_EB_xcorr_time_corrSac_onset( range_trials);
BEHAVE_EB_xcorr_time_corrSac_vmax   = BEHAVE_EB_xcorr_time_corrSac_vmax(  range_trials);
BEHAVE_EB_xcorr_time_corrSac_offset = BEHAVE_EB_xcorr_time_corrSac_offset(range_trials);

BEHAVE_EB_xcorr_time_cue_present(end+1)    = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_cue_present(end)])+1;
BEHAVE_EB_xcorr_time_primSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_primSac_onset(end)])+1;
BEHAVE_EB_xcorr_time_primSac_vmax(end+1)   = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_primSac_vmax(end)])+1;
BEHAVE_EB_xcorr_time_primSac_offset(end+1) = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_primSac_offset(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_onset(end+1)  = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_corrSac_onset(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_vmax(end+1)   = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_corrSac_vmax(end)])+1;
BEHAVE_EB_xcorr_time_corrSac_offset(end+1) = max([BEHAVE_EB_xcorr_time_1K(end), BEHAVE_EB_xcorr_time_corrSac_offset(end)])+1;
BEHAVE_EB_xcorr_ind_cue_present    = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_vmax   = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_primSac_offset = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_onset  = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_vmax   = nan(num_trials, 1);
BEHAVE_EB_xcorr_ind_corrSac_offset = nan(num_trials, 1);
counter_cue_present    = find(BEHAVE_EB_xcorr_time_cue_present    >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_primSac_onset  = find(BEHAVE_EB_xcorr_time_primSac_onset  >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_primSac_vmax   = find(BEHAVE_EB_xcorr_time_primSac_vmax   >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_primSac_offset = find(BEHAVE_EB_xcorr_time_primSac_offset >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_corrSac_onset  = find(BEHAVE_EB_xcorr_time_corrSac_onset  >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_corrSac_vmax   = find(BEHAVE_EB_xcorr_time_corrSac_vmax   >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_corrSac_offset = find(BEHAVE_EB_xcorr_time_corrSac_offset >= BEHAVE_EB_xcorr_time_1K(1), 1, 'first');
counter_trial_cue_present    = 1;
counter_trial_primSac_onset  = 1;
counter_trial_primSac_vmax   = 1;
counter_trial_primSac_offset = 1;
counter_trial_corrSac_onset  = 1;
counter_trial_corrSac_vmax   = 1;
counter_trial_corrSac_offset = 1;
for counter_time_point = 1 : length_time_
    time_ponit_ = BEHAVE_EB_xcorr_time_1K(counter_time_point);
    if time_ponit_>=BEHAVE_EB_xcorr_time_cue_present(counter_cue_present)
        BEHAVE_EB_xcorr_ind_cue_present(counter_trial_cue_present) = counter_time_point;
        counter_cue_present = counter_cue_present + 1;
        counter_trial_cue_present = counter_trial_cue_present + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_onset(counter_primSac_onset)
        BEHAVE_EB_xcorr_ind_primSac_onset(counter_trial_primSac_onset) = counter_time_point;
        counter_primSac_onset = counter_primSac_onset + 1;
        counter_trial_primSac_onset = counter_trial_primSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_vmax(counter_primSac_vmax)
        BEHAVE_EB_xcorr_ind_primSac_vmax(counter_trial_primSac_vmax) = counter_time_point;
        counter_primSac_vmax = counter_primSac_vmax + 1;
        counter_trial_primSac_vmax = counter_trial_primSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_primSac_offset(counter_primSac_offset)
        BEHAVE_EB_xcorr_ind_primSac_offset(counter_trial_primSac_offset) = counter_time_point;
        counter_primSac_offset = counter_primSac_offset + 1;
        counter_trial_primSac_offset = counter_trial_primSac_offset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_onset(counter_corrSac_onset)
        BEHAVE_EB_xcorr_ind_corrSac_onset(counter_trial_corrSac_onset) = counter_time_point;
        counter_corrSac_onset = counter_corrSac_onset + 1;
        counter_trial_corrSac_onset = counter_trial_corrSac_onset + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_vmax(counter_corrSac_vmax)
        BEHAVE_EB_xcorr_ind_corrSac_vmax(counter_trial_corrSac_vmax) = counter_time_point;
        counter_corrSac_vmax = counter_corrSac_vmax + 1;
        counter_trial_corrSac_vmax = counter_trial_corrSac_vmax + 1;
    end
    if time_ponit_>=BEHAVE_EB_xcorr_time_corrSac_offset(counter_corrSac_offset)
        BEHAVE_EB_xcorr_ind_corrSac_offset(counter_trial_corrSac_offset) = counter_time_point;
        counter_corrSac_offset = counter_corrSac_offset + 1;
        counter_trial_corrSac_offset = counter_trial_corrSac_offset + 1;
    end
end
% convert xcorr to aligned for EPHYS. We find the events on BEHAVE and then
% should convert it to EPHYS for SS & CS (EPHYS related events). We should not convert BEHAVE for
% BEHAVE related events.

% EPHYS_EB_aligned_ind_cue_present_1K    = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);
% EPHYS_EB_aligned_ind_cue_photodiode_1K = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);
EPHYS_EB_aligned_ind_cue_present_1K = EPHYS.CH_EVE.align_photodiode.EPHYS_PD_aligned_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);

EPHYS_EB_aligned_ind_primSac_onset_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_primSac_onset);
EPHYS_EB_aligned_ind_primSac_vmax_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_primSac_vmax);
EPHYS_EB_aligned_ind_primSac_offset_1K = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_primSac_offset);
EPHYS_EB_aligned_ind_corrSac_onset_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_onset);
EPHYS_EB_aligned_ind_corrSac_vmax_1K  = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_vmax);
EPHYS_EB_aligned_ind_corrSac_offset_1K = EPHYS.CH_EVE.align_states.EPHYS_EB_aligned_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_offset);
inds_span_cue_present    = ((-50+1) : 1 : (100))';
inds_span_primSac_onset  = ((-50+1) : 1 : (100))';
inds_span_primSac_vmax   = ((-50+1) : 1 : (100))';
inds_span_primSac_offset = ((-50+1) : 1 : (100))';
inds_span_corrSac_onset  = ((-50+1) : 1 : (100))';
inds_span_corrSac_vmax   = ((-50+1) : 1 : (100))';
inds_span_corrSac_offset = ((-50+1) : 1 : (100))';
% Build EPHYS_EB_aligned_inds
EPHYS_time_1K     = EPHYS.CH_EVE.EPHYS_time_1K;
length_time_      = length(EPHYS_time_1K);

EPHYS_EB_aligned_inds_cue_present_1K = repmat( EPHYS_EB_aligned_ind_cue_present_1K(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
EPHYS_EB_aligned_inds_cue_present_1K( EPHYS_EB_aligned_inds_cue_present_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_cue_present_1K( EPHYS_EB_aligned_inds_cue_present_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_onset_1K = repmat( EPHYS_EB_aligned_ind_primSac_onset_1K(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
EPHYS_EB_aligned_inds_primSac_onset_1K( EPHYS_EB_aligned_inds_primSac_onset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_onset_1K( EPHYS_EB_aligned_inds_primSac_onset_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_vmax_1K = repmat( EPHYS_EB_aligned_ind_primSac_vmax_1K(:), 1, length(inds_span_primSac_vmax)) + repmat(inds_span_primSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_primSac_vmax), 1);
EPHYS_EB_aligned_inds_primSac_vmax_1K( EPHYS_EB_aligned_inds_primSac_vmax_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_vmax_1K( EPHYS_EB_aligned_inds_primSac_vmax_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_primSac_offset_1K = repmat( EPHYS_EB_aligned_ind_primSac_offset_1K(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
EPHYS_EB_aligned_inds_primSac_offset_1K( EPHYS_EB_aligned_inds_primSac_offset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_primSac_offset_1K( EPHYS_EB_aligned_inds_primSac_offset_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_onset_1K = repmat( EPHYS_EB_aligned_ind_corrSac_onset_1K(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
EPHYS_EB_aligned_inds_corrSac_onset_1K( EPHYS_EB_aligned_inds_corrSac_onset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_onset_1K( EPHYS_EB_aligned_inds_corrSac_onset_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_vmax_1K = repmat( EPHYS_EB_aligned_ind_corrSac_vmax_1K(:), 1, length(inds_span_corrSac_vmax)) + repmat(inds_span_corrSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_corrSac_vmax), 1);
EPHYS_EB_aligned_inds_corrSac_vmax_1K( EPHYS_EB_aligned_inds_corrSac_vmax_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_vmax_1K( EPHYS_EB_aligned_inds_corrSac_vmax_1K > length_time_ ) = length_time_;

EPHYS_EB_aligned_inds_corrSac_offset_1K = repmat( EPHYS_EB_aligned_ind_corrSac_offset_1K(:), 1, length(inds_span_corrSac_offset)) + repmat(inds_span_corrSac_offset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_offset), 1);
EPHYS_EB_aligned_inds_corrSac_offset_1K( EPHYS_EB_aligned_inds_corrSac_offset_1K < 1 ) = 1;
EPHYS_EB_aligned_inds_corrSac_offset_1K( EPHYS_EB_aligned_inds_corrSac_offset_1K > length_time_ ) = length_time_;

% We should not convert BEHAVE for BEHAVE related events.
BEHAVE_EB_aligned_ind_cue_present_1K    = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_cue_present);
BEHAVE_EB_aligned_ind_primSac_onset_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_primSac_onset);
BEHAVE_EB_aligned_ind_primSac_vmax_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_primSac_vmax);
BEHAVE_EB_aligned_ind_primSac_offset_1K = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_primSac_offset);
BEHAVE_EB_aligned_ind_corrSac_onset_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_onset);
BEHAVE_EB_aligned_ind_corrSac_vmax_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_vmax);
BEHAVE_EB_aligned_ind_corrSac_offset_1K  = EPHYS.CH_EVE.align_states.BEHAVE_EB_xcorr_ind_1K(BEHAVE_EB_xcorr_ind_corrSac_offset);
% Build BEHAVE_EB_aligned_inds
BEHAVE_time_1K     = EPHYS.CH_EVE.BEHAVE_time_1K;
length_time_      = length(BEHAVE_time_1K);

BEHAVE_EB_aligned_inds_cue_present_1K = repmat( BEHAVE_EB_aligned_ind_cue_present_1K(:), 1, length(inds_span_cue_present)) + repmat(inds_span_cue_present(:)', length(BEHAVE_EB_xcorr_ind_cue_present), 1);
BEHAVE_EB_aligned_inds_cue_present_1K( BEHAVE_EB_aligned_inds_cue_present_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_cue_present_1K( BEHAVE_EB_aligned_inds_cue_present_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_onset_1K = repmat( BEHAVE_EB_aligned_ind_primSac_onset_1K(:), 1, length(inds_span_primSac_onset)) + repmat(inds_span_primSac_onset(:)', length(BEHAVE_EB_xcorr_ind_primSac_onset), 1);
BEHAVE_EB_aligned_inds_primSac_onset_1K( BEHAVE_EB_aligned_inds_primSac_onset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_onset_1K( BEHAVE_EB_aligned_inds_primSac_onset_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_vmax_1K = repmat( BEHAVE_EB_aligned_ind_primSac_vmax_1K(:), 1, length(inds_span_primSac_vmax)) + repmat(inds_span_primSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_primSac_vmax), 1);
BEHAVE_EB_aligned_inds_primSac_vmax_1K( BEHAVE_EB_aligned_inds_primSac_vmax_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_vmax_1K( BEHAVE_EB_aligned_inds_primSac_vmax_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_primSac_offset_1K = repmat( BEHAVE_EB_aligned_ind_primSac_offset_1K(:), 1, length(inds_span_primSac_offset)) + repmat(inds_span_primSac_offset(:)', length(BEHAVE_EB_xcorr_ind_primSac_offset), 1);
BEHAVE_EB_aligned_inds_primSac_offset_1K( BEHAVE_EB_aligned_inds_primSac_offset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_primSac_offset_1K( BEHAVE_EB_aligned_inds_primSac_offset_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_onset_1K = repmat( BEHAVE_EB_aligned_ind_corrSac_onset_1K(:), 1, length(inds_span_corrSac_onset)) + repmat(inds_span_corrSac_onset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_onset), 1);
BEHAVE_EB_aligned_inds_corrSac_onset_1K( BEHAVE_EB_aligned_inds_corrSac_onset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_onset_1K( BEHAVE_EB_aligned_inds_corrSac_onset_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_vmax_1K = repmat( BEHAVE_EB_aligned_ind_corrSac_vmax_1K(:), 1, length(inds_span_corrSac_vmax)) + repmat(inds_span_corrSac_vmax(:)', length(BEHAVE_EB_xcorr_ind_corrSac_vmax), 1);
BEHAVE_EB_aligned_inds_corrSac_vmax_1K( BEHAVE_EB_aligned_inds_corrSac_vmax_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_vmax_1K( BEHAVE_EB_aligned_inds_corrSac_vmax_1K > length_time_ ) = length_time_;

BEHAVE_EB_aligned_inds_corrSac_offset_1K = repmat( BEHAVE_EB_aligned_ind_corrSac_offset_1K(:), 1, length(inds_span_corrSac_offset)) + repmat(inds_span_corrSac_offset(:)', length(BEHAVE_EB_xcorr_ind_corrSac_offset), 1);
BEHAVE_EB_aligned_inds_corrSac_offset_1K( BEHAVE_EB_aligned_inds_corrSac_offset_1K < 1 ) = 1;
BEHAVE_EB_aligned_inds_corrSac_offset_1K( BEHAVE_EB_aligned_inds_corrSac_offset_1K > length_time_ ) = length_time_;

EPHYS.CH_EVE.EPHYS_EB_aligned_ind_cue_present_1K     = EPHYS_EB_aligned_ind_cue_present_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_onset_1K   = EPHYS_EB_aligned_ind_primSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_vmax_1K    = EPHYS_EB_aligned_ind_primSac_vmax_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_offset_1K  = EPHYS_EB_aligned_ind_primSac_offset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_onset_1K   = EPHYS_EB_aligned_ind_corrSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_vmax_1K    = EPHYS_EB_aligned_ind_corrSac_vmax_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_offset_1K  = EPHYS_EB_aligned_ind_corrSac_offset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_cue_present_1K    = EPHYS_EB_aligned_inds_cue_present_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_onset_1K  = EPHYS_EB_aligned_inds_primSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_vmax_1K   = EPHYS_EB_aligned_inds_primSac_vmax_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_primSac_offset_1K = EPHYS_EB_aligned_inds_primSac_offset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_onset_1K  = EPHYS_EB_aligned_inds_corrSac_onset_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_vmax_1K   = EPHYS_EB_aligned_inds_corrSac_vmax_1K;
EPHYS.CH_EVE.EPHYS_EB_aligned_inds_corrSac_offset_1K = EPHYS_EB_aligned_inds_corrSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_cue_present_1K     = BEHAVE_EB_aligned_ind_cue_present_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_onset_1K   = BEHAVE_EB_aligned_ind_primSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_vmax_1K    = BEHAVE_EB_aligned_ind_primSac_vmax_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_primSac_offset_1K  = BEHAVE_EB_aligned_ind_primSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_onset_1K   = BEHAVE_EB_aligned_ind_corrSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_vmax_1K    = BEHAVE_EB_aligned_ind_corrSac_vmax_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_ind_corrSac_offset_1K  = BEHAVE_EB_aligned_ind_corrSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_cue_present_1K    = BEHAVE_EB_aligned_inds_cue_present_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_onset_1K  = BEHAVE_EB_aligned_inds_primSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_vmax_1K   = BEHAVE_EB_aligned_inds_primSac_vmax_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_primSac_offset_1K = BEHAVE_EB_aligned_inds_primSac_offset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_onset_1K  = BEHAVE_EB_aligned_inds_corrSac_onset_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_vmax_1K   = BEHAVE_EB_aligned_inds_corrSac_vmax_1K;
EPHYS.CH_EVE.BEHAVE_EB_aligned_inds_corrSac_offset_1K = BEHAVE_EB_aligned_inds_corrSac_offset_1K;
EPHYS.CH_EVE.inds_span_cue_present      = inds_span_cue_present(:)';
EPHYS.CH_EVE.inds_span_primSac_onset    = inds_span_primSac_onset(:)';
EPHYS.CH_EVE.inds_span_primSac_vmax     = inds_span_primSac_vmax(:)';
EPHYS.CH_EVE.inds_span_primSac_offset   = inds_span_primSac_offset(:)';
EPHYS.CH_EVE.inds_span_corrSac_onset    = inds_span_corrSac_onset(:)';
EPHYS.CH_EVE.inds_span_corrSac_vmax     = inds_span_corrSac_vmax(:)';
EPHYS.CH_EVE.inds_span_corrSac_offset   = inds_span_corrSac_offset(:)';

fprintf(' --> Completed. \n')

end

%% function concatenate_dataset
function upper_field_struct = concatenate_dataset(dataset_, upper_field_name, horz_OR_vert)
if ~isempty(upper_field_name)
dataset = struct(upper_field_name,struct());
field_names_ = fieldnames(dataset_(1).(upper_field_name));
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(upper_field_name).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset.(upper_field_name), field_names_{counter_fields})
            dataset.(upper_field_name).(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(upper_field_name).(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(upper_field_name).(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset.(upper_field_name);
elseif isempty(upper_field_name)
dataset = struct();
field_names_ = fieldnames(dataset_);
for counter_fields = 1 : 1 : length(field_names_)
    for counter_dataset = 1 : 1 : length(dataset_)
        variable_TRIALS_DATA_ALL_ = dataset_(counter_dataset).(field_names_{counter_fields});
        % the field does not exist in TRIALS_DATA
        if ~isfield(dataset, field_names_{counter_fields})
            dataset.(field_names_{counter_fields}) = [];
        end
        variable_TRIALS_DATA_ = dataset.(field_names_{counter_fields});
        variable_TRIALS_DATA_ = horz_OR_vert(variable_TRIALS_DATA_, variable_TRIALS_DATA_ALL_);
        dataset.(field_names_{counter_fields}) = variable_TRIALS_DATA_;
    end
end
upper_field_struct = dataset;
end
end

%% function plot_rasters_data
function fig_handle_ = plot_rasters_data(raster_data, plot_data)
%%
fig_num_               = plot_data.fig_num_;
xlabel_text_raster_    = plot_data.xlabel_text_raster_;
xlabel_text_CS_probab_ = plot_data.xlabel_text_CS_probab_;
inds_span              = plot_data.inds_span;

train_data_logic_SS_000 = raster_data.train_data_logic_SS_000;
train_data_logic_CS_000 = raster_data.train_data_logic_CS_000;
velocity_data_000 = raster_data.velocity_data_000;
train_data_logic_cue_present_000    = raster_data.train_data_logic_cue_present_000;
train_data_logic_primSac_onset_000  = raster_data.train_data_logic_primSac_onset_000;
train_data_logic_primSac_vmax_000   = raster_data.train_data_logic_primSac_vmax_000;
train_data_logic_primSac_offset_000 = raster_data.train_data_logic_primSac_offset_000;
train_data_logic_corrSac_onset_000  = raster_data.train_data_logic_corrSac_onset_000;
train_data_logic_corrSac_vmax_000   = raster_data.train_data_logic_corrSac_vmax_000;
train_data_logic_corrSac_offset_000 = raster_data.train_data_logic_corrSac_offset_000;

train_data_logic_SS_045 = raster_data.train_data_logic_SS_045;
train_data_logic_CS_045 = raster_data.train_data_logic_CS_045;
velocity_data_045 = raster_data.velocity_data_045;
train_data_logic_cue_present_045    = raster_data.train_data_logic_cue_present_045;
train_data_logic_primSac_onset_045  = raster_data.train_data_logic_primSac_onset_045;
train_data_logic_primSac_vmax_045   = raster_data.train_data_logic_primSac_vmax_045;
train_data_logic_primSac_offset_045 = raster_data.train_data_logic_primSac_offset_045;
train_data_logic_corrSac_onset_045  = raster_data.train_data_logic_corrSac_onset_045;
train_data_logic_corrSac_vmax_045   = raster_data.train_data_logic_corrSac_vmax_045;
train_data_logic_corrSac_offset_045 = raster_data.train_data_logic_corrSac_offset_045;

train_data_logic_SS_090 = raster_data.train_data_logic_SS_090;
train_data_logic_CS_090 = raster_data.train_data_logic_CS_090;
velocity_data_090 = raster_data.velocity_data_090;
train_data_logic_cue_present_090    = raster_data.train_data_logic_cue_present_090;
train_data_logic_primSac_onset_090  = raster_data.train_data_logic_primSac_onset_090;
train_data_logic_primSac_vmax_090   = raster_data.train_data_logic_primSac_vmax_090;
train_data_logic_primSac_offset_090 = raster_data.train_data_logic_primSac_offset_090;
train_data_logic_corrSac_onset_090  = raster_data.train_data_logic_corrSac_onset_090;
train_data_logic_corrSac_vmax_090   = raster_data.train_data_logic_corrSac_vmax_090;
train_data_logic_corrSac_offset_090 = raster_data.train_data_logic_corrSac_offset_090;

train_data_logic_SS_135 = raster_data.train_data_logic_SS_135;
train_data_logic_CS_135 = raster_data.train_data_logic_CS_135;
velocity_data_135 = raster_data.velocity_data_135;
train_data_logic_cue_present_135    = raster_data.train_data_logic_cue_present_135;
train_data_logic_primSac_onset_135  = raster_data.train_data_logic_primSac_onset_135;
train_data_logic_primSac_vmax_135   = raster_data.train_data_logic_primSac_vmax_135;
train_data_logic_primSac_offset_135 = raster_data.train_data_logic_primSac_offset_135;
train_data_logic_corrSac_onset_135  = raster_data.train_data_logic_corrSac_onset_135;
train_data_logic_corrSac_vmax_135   = raster_data.train_data_logic_corrSac_vmax_135;
train_data_logic_corrSac_offset_135 = raster_data.train_data_logic_corrSac_offset_135;

train_data_logic_SS_180 = raster_data.train_data_logic_SS_180;
train_data_logic_CS_180 = raster_data.train_data_logic_CS_180;
velocity_data_180 = raster_data.velocity_data_180;
train_data_logic_cue_present_180    = raster_data.train_data_logic_cue_present_180;
train_data_logic_primSac_onset_180  = raster_data.train_data_logic_primSac_onset_180;
train_data_logic_primSac_vmax_180   = raster_data.train_data_logic_primSac_vmax_180;
train_data_logic_primSac_offset_180 = raster_data.train_data_logic_primSac_offset_180;
train_data_logic_corrSac_onset_180  = raster_data.train_data_logic_corrSac_onset_180;
train_data_logic_corrSac_vmax_180   = raster_data.train_data_logic_corrSac_vmax_180;
train_data_logic_corrSac_offset_180 = raster_data.train_data_logic_corrSac_offset_180;

train_data_logic_SS_225 = raster_data.train_data_logic_SS_225;
train_data_logic_CS_225 = raster_data.train_data_logic_CS_225;
velocity_data_225 = raster_data.velocity_data_225;
train_data_logic_cue_present_225    = raster_data.train_data_logic_cue_present_225;
train_data_logic_primSac_onset_225  = raster_data.train_data_logic_primSac_onset_225;
train_data_logic_primSac_vmax_225   = raster_data.train_data_logic_primSac_vmax_225;
train_data_logic_primSac_offset_225 = raster_data.train_data_logic_primSac_offset_225;
train_data_logic_corrSac_onset_225  = raster_data.train_data_logic_corrSac_onset_225;
train_data_logic_corrSac_vmax_225   = raster_data.train_data_logic_corrSac_vmax_225;
train_data_logic_corrSac_offset_225 = raster_data.train_data_logic_corrSac_offset_225;

train_data_logic_SS_270 = raster_data.train_data_logic_SS_270;
train_data_logic_CS_270 = raster_data.train_data_logic_CS_270;
velocity_data_270 = raster_data.velocity_data_270;
train_data_logic_cue_present_270    = raster_data.train_data_logic_cue_present_270;
train_data_logic_primSac_onset_270  = raster_data.train_data_logic_primSac_onset_270;
train_data_logic_primSac_vmax_270   = raster_data.train_data_logic_primSac_vmax_270;
train_data_logic_primSac_offset_270 = raster_data.train_data_logic_primSac_offset_270;
train_data_logic_corrSac_onset_270  = raster_data.train_data_logic_corrSac_onset_270;
train_data_logic_corrSac_vmax_270   = raster_data.train_data_logic_corrSac_vmax_270;
train_data_logic_corrSac_offset_270 = raster_data.train_data_logic_corrSac_offset_270;

train_data_logic_SS_315 = raster_data.train_data_logic_SS_315;
train_data_logic_CS_315 = raster_data.train_data_logic_CS_315;
velocity_data_315 = raster_data.velocity_data_315;
train_data_logic_cue_present_315    = raster_data.train_data_logic_cue_present_315;
train_data_logic_primSac_onset_315  = raster_data.train_data_logic_primSac_onset_315;
train_data_logic_primSac_vmax_315   = raster_data.train_data_logic_primSac_vmax_315;
train_data_logic_primSac_offset_315 = raster_data.train_data_logic_primSac_offset_315;
train_data_logic_corrSac_onset_315  = raster_data.train_data_logic_corrSac_onset_315;
train_data_logic_corrSac_vmax_315   = raster_data.train_data_logic_corrSac_vmax_315;
train_data_logic_corrSac_offset_315 = raster_data.train_data_logic_corrSac_offset_315;

% % CS Probab
% if contains(xlabel_text_CS_probab_, '+200')
%     range_inds_probability = 301:500;
% elseif contains(xlabel_text_CS_probab_, '-200')
%     range_inds_probability = 101:300;
% end
% prob_000  = sum( sum(train_data_logic_SS_000(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_000, 1);
% prob_045  = sum( sum(train_data_logic_SS_045(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_045, 1);
% prob_090  = sum( sum(train_data_logic_SS_090(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_090, 1);
% prob_135  = sum( sum(train_data_logic_SS_135(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_135, 1);
% prob_180  = sum( sum(train_data_logic_SS_180(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_180, 1);
% prob_225  = sum( sum(train_data_logic_SS_225(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_225, 1);
% prob_270  = sum( sum(train_data_logic_SS_270(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_270, 1);
% prob_315  = sum( sum(train_data_logic_SS_315(:,range_inds_probability) ,2) > 0 ) / size(train_data_logic_SS_315, 1);
% prob_amplitude = [prob_000 prob_045 prob_090 prob_135 prob_180 prob_225 prob_270 prob_315 prob_000];
% % plot xlim and ylim
range_SS_Firing = [0 200];
range_vm        = [0 600];
% % xlabel text
% xlabel_text_raster_ = {'Time relative to cue presentation (ms)', 'Directions based on primary sac'};
% xlabel_text_CS_probab_ = {'CS probability based on [0 +200]ms'};
% % figure
% fig_num_ = 1;
Line_Color = lines(7);
color_SS = Line_Color(1,:);
color_SS_firing = [0    0.3    0.5];
color_CS = Line_Color(7,:);
color_cue_present = Line_Color(3,:);
color_primSac_onset = Line_Color(4,:);
color_primSac_vmax = Line_Color(5,:);
color_primSac_offset = Line_Color(2,:);
color_corrSac_onset = Line_Color(4,:);
color_corrSac_vmax = Line_Color(5,:);
color_corrSac_offset = Line_Color(2,:);

fig_handle_ = figure(fig_num_);
clf(fig_handle_)

%% 000
subplot(3,3,6)
hold on
train_data_logic_SS_ = train_data_logic_SS_000;
train_data_logic_CS_ = train_data_logic_CS_000;
train_data_logic_cue_present_    = train_data_logic_cue_present_000;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_000;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_000;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_000;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_000;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_000;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_000;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)

xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_000), '-g', 'LineWidth' , 2)
ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 045
subplot(3,3,3)
hold on
train_data_logic_SS_ = train_data_logic_SS_045;
train_data_logic_CS_ = train_data_logic_CS_045;
train_data_logic_cue_present_    = train_data_logic_cue_present_045;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_045;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_045;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_045;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_045;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_045;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_045;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_045), '-g', 'LineWidth' , 2)
ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 090
subplot(3,3,2)
hold on
train_data_logic_SS_ = train_data_logic_SS_090;
train_data_logic_CS_ = train_data_logic_CS_090;
train_data_logic_cue_present_    = train_data_logic_cue_present_090;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_090;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_090;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_090;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_090;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_090;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_090;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_090), '-g', 'LineWidth' , 2)
% ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 135
subplot(3,3,1)
hold on
train_data_logic_SS_ = train_data_logic_SS_135;
train_data_logic_CS_ = train_data_logic_CS_135;
train_data_logic_cue_present_    = train_data_logic_cue_present_135;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_135;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_135;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_135;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_135;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_135;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_135;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_135), '-g', 'LineWidth' , 2)
% ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 180
subplot(3,3,4)
hold on
train_data_logic_SS_ = train_data_logic_SS_180;
train_data_logic_CS_ = train_data_logic_CS_180;
train_data_logic_cue_present_    = train_data_logic_cue_present_180;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_180;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_180;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_180;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_180;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_180;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_180;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_180), '-g', 'LineWidth' , 2)
% ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 225
subplot(3,3,7)
hold on
train_data_logic_SS_ = train_data_logic_SS_225;
train_data_logic_CS_ = train_data_logic_CS_225;
train_data_logic_cue_present_    = train_data_logic_cue_present_225;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_225;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_225;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_225;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_225;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_225;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_225;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_225), '-g', 'LineWidth' , 2)
% ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 270
subplot(3,3,8)
hold on
train_data_logic_SS_ = train_data_logic_SS_270;
train_data_logic_CS_ = train_data_logic_CS_270;
train_data_logic_cue_present_    = train_data_logic_cue_present_270;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_270;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_270;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_270;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_270;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_270;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_270;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_270), '-g', 'LineWidth' , 2)
% ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% 315
subplot(3,3,9)
hold on
train_data_logic_SS_ = train_data_logic_SS_315;
train_data_logic_CS_ = train_data_logic_CS_315;
train_data_logic_cue_present_    = train_data_logic_cue_present_315;
train_data_logic_primSac_onset_  = train_data_logic_primSac_onset_315;
train_data_logic_primSac_vmax_   = train_data_logic_primSac_vmax_315;
train_data_logic_primSac_offset_ = train_data_logic_primSac_offset_315;
train_data_logic_corrSac_onset_  = train_data_logic_corrSac_onset_315;
train_data_logic_corrSac_vmax_   = train_data_logic_corrSac_vmax_315;
train_data_logic_corrSac_offset_ = train_data_logic_corrSac_offset_315;

[x_axis_SS_, y_axis_SS_] = ESN_raster_plot_axes(train_data_logic_SS_, inds_span, 0.5);
plot(x_axis_SS_(:), y_axis_SS_(:), 'LineWidth', 1, 'Color', color_SS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_CS_, inds_span, 1);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 3, 'Color', color_CS)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_cue_present_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_cue_present)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_onset)
% [x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_vmax_, inds_span, 0.5);
% plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_primSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_primSac_offset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_onset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_onset)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_vmax_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_vmax)
[x_axis_CS_, y_axis_CS_] = ESN_raster_plot_axes(train_data_logic_corrSac_offset_, inds_span, 0.5);
plot(x_axis_CS_(:), y_axis_CS_(:), 'LineWidth', 2, 'Color', color_corrSac_offset)
xlim([min(inds_span)-1 max(inds_span)+1])
ylim([(1-3) (size(train_data_logic_CS_,1)+3)])
% ylabel('Trials')

yyaxis right;
firing_SS_ = mean(train_data_logic_SS_) * 1000;
plot(inds_span, ESN_smooth(firing_SS_), 'LineWidth', 2, 'Color', color_SS_firing)
plot(inds_span, nanmean(velocity_data_315), '-g', 'LineWidth' , 2)
ylabel('SS Firing (spk/s)')
%ylim(range_SS_Firing)
ylim([0 50])
set(gca, 'YColor', color_SS_firing)
xlim([min(inds_span)-1 max(inds_span)+1])

%% Probability
% subplot(3,3,5)
% hold on
% 
% plot(cosd(0:5:360)*0.10, sind(0:5:360)*0.10, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% plot(cosd(0:5:360)*0.20, sind(0:5:360)*0.20, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% plot(cosd(0:5:360)*0.30, sind(0:5:360)*0.30, 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
% 
% x_axis = [cosd(0) cosd(45) cosd(90) cosd(135) cosd(180) cosd(225) cosd(270) cosd(315) cosd(0)] .* prob_amplitude;
% y_axis = [sind(0) sind(45) sind(90) sind(135) sind(180) sind(225) sind(270) sind(315) sind(0)] .* prob_amplitude;
% x_axis = x_axis(~isnan(x_axis)); y_axis = y_axis(~isnan(y_axis));
% plot(x_axis(:), y_axis(:), 'LineWidth', 3, 'Color', color_CS)
% axis equal
% set(gca, 'YColor', color_CS)
% set(gca, 'XColor', color_CS)

%% xlabel
subplot(3,3,5)
xlabel(xlabel_text_CS_probab_);

subplot(3,3,8)
xlabel(xlabel_text_raster_);

 ESN_Beautify_Plot(fig_handle_, [8.0 8.0])

end

%% function single_dataset_raster
function raster_data = single_dataset_raster(EPHYS, BEHAVE, EPHYS_inds_event, BEHAVE_inds_event, prim_OR_corr)

% % CS and SS data
SS_train_aligned     = EPHYS.CH_EVE.EPHYS_SS_train_1K;
CS_train_aligned     = EPHYS.CH_EVE.EPHYS_CS_train_1K;

EPHYS_cue_present_train_aligned    = false(size(SS_train_aligned));
EPHYS_primSac_onset_train_aligned  = false(size(SS_train_aligned));
EPHYS_primSac_vmax_train_aligned   = false(size(SS_train_aligned));
EPHYS_primSac_offset_train_aligned = false(size(SS_train_aligned));
EPHYS_corrSac_onset_train_aligned  = false(size(SS_train_aligned));
EPHYS_corrSac_vmax_train_aligned   = false(size(SS_train_aligned));
EPHYS_corrSac_offset_train_aligned = false(size(SS_train_aligned));

EPHYS_cue_present_train_aligned( EPHYS.CH_EVE.EPHYS_EB_aligned_ind_cue_present_1K)      = true;
EPHYS_primSac_onset_train_aligned( EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_onset_1K)  = true;
EPHYS_primSac_vmax_train_aligned(  EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_vmax_1K)   = true;
EPHYS_primSac_offset_train_aligned(EPHYS.CH_EVE.EPHYS_EB_aligned_ind_primSac_offset_1K) = true;
EPHYS_corrSac_onset_train_aligned( EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_onset_1K)  = true;
EPHYS_corrSac_vmax_train_aligned(  EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_vmax_1K)   = true;
EPHYS_corrSac_offset_train_aligned(EPHYS.CH_EVE.EPHYS_EB_aligned_ind_corrSac_offset_1K) = true;

% % eye velocity
BEHAVE_eye_r_vm_filt = EPHYS.CH_EVE.BEHAVE_eye_r_vm_filt_1K;
% % check validity of trials
num_trials = size(BEHAVE_inds_event,1);
% inds_valid = BEHAVE.SACS_PRIM_DATA.validity(1:num_trials) & BEHAVE.SACS_CORR_DATA.validity(1:num_trials);
% error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish(1:num_trials) - BEHAVE.TRIALS_DATA.cue_x(1:num_trials)).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish(1:num_trials) - BEHAVE.TRIALS_DATA.cue_y(1:num_trials)).^2 );
% error_corr = sqrt( (BEHAVE.SACS_CORR_DATA.eye_r_px_finish(1:num_trials) - BEHAVE.TRIALS_DATA.end_x(1:num_trials)).^2 + (BEHAVE.SACS_CORR_DATA.eye_r_py_finish(1:num_trials) - BEHAVE.TRIALS_DATA.end_y(1:num_trials)).^2 );
% inds_valid = inds_valid & (error_prim<3) & (error_corr<3);
inds_valid = BEHAVE.SACS_PRIM_DATA.validity(1:num_trials);
error_prim = sqrt( (BEHAVE.SACS_PRIM_DATA.eye_r_px_finish(1:num_trials) - BEHAVE.TRIALS_DATA.cue_x(1:num_trials)).^2 + (BEHAVE.SACS_PRIM_DATA.eye_r_py_finish(1:num_trials) - BEHAVE.TRIALS_DATA.cue_y(1:num_trials)).^2 );
inds_valid = inds_valid & (error_prim<3);
% % inds directions
if contains(prim_OR_corr, 'prim')
prim_delta_x = (BEHAVE.TRIALS_DATA.cue_x(1:num_trials)-BEHAVE.TRIALS_DATA.start_x(1:num_trials));
prim_delta_y = (BEHAVE.TRIALS_DATA.cue_y(1:num_trials)-BEHAVE.TRIALS_DATA.start_y(1:num_trials));
prim_angle = atan2d(prim_delta_y, prim_delta_x);
inds_000   = (prim_angle >  -20) & (prim_angle <  +20) & inds_valid;
inds_045   = (prim_angle >  +25) & (prim_angle <  +65) & inds_valid;
inds_090   = (prim_angle >  +70) & (prim_angle < +110) & inds_valid;
inds_135   = (prim_angle > +115) & (prim_angle < +155) & inds_valid;
inds_225   = (prim_angle > -155) & (prim_angle < -115) & inds_valid;
inds_270   = (prim_angle > -110) & (prim_angle <  -70) & inds_valid;
inds_315   = (prim_angle >  -65) & (prim_angle <  -25) & inds_valid;
inds_180   = ((prim_angle > +155) | (prim_angle < -155)) & inds_valid;
elseif contains(prim_OR_corr, 'corr')
corr_delta_x = (BEHAVE.TRIALS_DATA.end_x(1:num_trials)-BEHAVE.TRIALS_DATA.cue_x(1:num_trials));
corr_delta_y = (BEHAVE.TRIALS_DATA.end_y(1:num_trials)-BEHAVE.TRIALS_DATA.cue_y(1:num_trials));
corr_angle = atan2d(corr_delta_y, corr_delta_x);
inds_000   = (corr_angle >  -20) & (corr_angle <  +20) & inds_valid;
inds_045   = (corr_angle >  +25) & (corr_angle <  +65) & inds_valid;
inds_090   = (corr_angle >  +70) & (corr_angle < +110) & inds_valid;
inds_135   = (corr_angle > +115) & (corr_angle < +155) & inds_valid;
inds_225   = (corr_angle > -155) & (corr_angle < -115) & inds_valid;
inds_270   = (corr_angle > -110) & (corr_angle <  -70) & inds_valid;
inds_315   = (corr_angle >  -65) & (corr_angle <  -25) & inds_valid;
inds_180   = ((corr_angle > +155) | (corr_angle < -155)) & inds_valid;
end

train_data_logic_SS_000 = SS_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_SS_045 = SS_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_SS_090 = SS_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_SS_135 = SS_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_SS_180 = SS_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_SS_225 = SS_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_SS_270 = SS_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_SS_315 = SS_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_CS_000 = CS_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_CS_045 = CS_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_CS_090 = CS_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_CS_135 = CS_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_CS_180 = CS_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_CS_225 = CS_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_CS_270 = CS_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_CS_315 = CS_train_aligned(EPHYS_inds_event(inds_315,:));

velocity_data_000 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_000,:));
velocity_data_045 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_045,:));
velocity_data_090 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_090,:));
velocity_data_135 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_135,:));
velocity_data_180 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_180,:));
velocity_data_225 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_225,:));
velocity_data_270 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_270,:));
velocity_data_315 = BEHAVE_eye_r_vm_filt(BEHAVE_inds_event(inds_315,:));

train_data_logic_cue_present_000 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_cue_present_045 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_cue_present_090 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_cue_present_135 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_cue_present_180 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_cue_present_225 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_cue_present_270 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_cue_present_315 = EPHYS_cue_present_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_primSac_onset_000 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_primSac_onset_045 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_primSac_onset_090 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_primSac_onset_135 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_primSac_onset_180 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_primSac_onset_225 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_primSac_onset_270 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_primSac_onset_315 = EPHYS_primSac_onset_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_primSac_vmax_000 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_primSac_vmax_045 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_primSac_vmax_090 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_primSac_vmax_135 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_primSac_vmax_180 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_primSac_vmax_225 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_primSac_vmax_270 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_primSac_vmax_315 = EPHYS_primSac_vmax_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_primSac_offset_000 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_primSac_offset_045 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_primSac_offset_090 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_primSac_offset_135 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_primSac_offset_180 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_primSac_offset_225 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_primSac_offset_270 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_primSac_offset_315 = EPHYS_primSac_offset_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_corrSac_onset_000 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_corrSac_onset_045 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_corrSac_onset_090 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_corrSac_onset_135 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_corrSac_onset_180 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_corrSac_onset_225 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_corrSac_onset_270 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_corrSac_onset_315 = EPHYS_corrSac_onset_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_corrSac_vmax_000 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_corrSac_vmax_045 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_corrSac_vmax_090 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_corrSac_vmax_135 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_corrSac_vmax_180 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_corrSac_vmax_225 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_corrSac_vmax_270 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_corrSac_vmax_315 = EPHYS_corrSac_vmax_train_aligned(EPHYS_inds_event(inds_315,:));

train_data_logic_corrSac_offset_000 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_000,:));
train_data_logic_corrSac_offset_045 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_045,:));
train_data_logic_corrSac_offset_090 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_090,:));
train_data_logic_corrSac_offset_135 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_135,:));
train_data_logic_corrSac_offset_180 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_180,:));
train_data_logic_corrSac_offset_225 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_225,:));
train_data_logic_corrSac_offset_270 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_270,:));
train_data_logic_corrSac_offset_315 = EPHYS_corrSac_offset_train_aligned(EPHYS_inds_event(inds_315,:));

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_000, 2)
raster_data.train_data_logic_SS_000 = train_data_logic_SS_000;
raster_data.train_data_logic_CS_000 = train_data_logic_CS_000;
raster_data.velocity_data_000       = velocity_data_000;
raster_data.train_data_logic_cue_present_000    = train_data_logic_cue_present_000;
raster_data.train_data_logic_primSac_onset_000  = train_data_logic_primSac_onset_000;
raster_data.train_data_logic_primSac_vmax_000   = train_data_logic_primSac_vmax_000;
raster_data.train_data_logic_primSac_offset_000 = train_data_logic_primSac_offset_000;
raster_data.train_data_logic_corrSac_onset_000  = train_data_logic_corrSac_onset_000;
raster_data.train_data_logic_corrSac_vmax_000   = train_data_logic_corrSac_vmax_000;
raster_data.train_data_logic_corrSac_offset_000 = train_data_logic_corrSac_offset_000;
else
raster_data.train_data_logic_SS_000 = train_data_logic_SS_000';
raster_data.train_data_logic_CS_000 = train_data_logic_CS_000';
raster_data.velocity_data_000       = velocity_data_000';
raster_data.train_data_logic_cue_present_000    = train_data_logic_cue_present_000';
raster_data.train_data_logic_primSac_onset_000  = train_data_logic_primSac_onset_000';
raster_data.train_data_logic_primSac_vmax_000   = train_data_logic_primSac_vmax_000';
raster_data.train_data_logic_primSac_offset_000 = train_data_logic_primSac_offset_000';
raster_data.train_data_logic_corrSac_onset_000  = train_data_logic_corrSac_onset_000';
raster_data.train_data_logic_corrSac_vmax_000   = train_data_logic_corrSac_vmax_000';
raster_data.train_data_logic_corrSac_offset_000 = train_data_logic_corrSac_offset_000';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_045, 2)
raster_data.train_data_logic_SS_045 = train_data_logic_SS_045;
raster_data.train_data_logic_CS_045 = train_data_logic_CS_045;
raster_data.velocity_data_045       = velocity_data_045;
raster_data.train_data_logic_cue_present_045    = train_data_logic_cue_present_045;
raster_data.train_data_logic_primSac_onset_045  = train_data_logic_primSac_onset_045;
raster_data.train_data_logic_primSac_vmax_045   = train_data_logic_primSac_vmax_045;
raster_data.train_data_logic_primSac_offset_045 = train_data_logic_primSac_offset_045;
raster_data.train_data_logic_corrSac_onset_045  = train_data_logic_corrSac_onset_045;
raster_data.train_data_logic_corrSac_vmax_045   = train_data_logic_corrSac_vmax_045;
raster_data.train_data_logic_corrSac_offset_045 = train_data_logic_corrSac_offset_045;
else
raster_data.train_data_logic_SS_045 = train_data_logic_SS_045';
raster_data.train_data_logic_CS_045 = train_data_logic_CS_045';
raster_data.velocity_data_045       = velocity_data_045';
raster_data.train_data_logic_cue_present_045    = train_data_logic_cue_present_045';
raster_data.train_data_logic_primSac_onset_045  = train_data_logic_primSac_onset_045';
raster_data.train_data_logic_primSac_vmax_045   = train_data_logic_primSac_vmax_045';
raster_data.train_data_logic_primSac_offset_045 = train_data_logic_primSac_offset_045';
raster_data.train_data_logic_corrSac_onset_045  = train_data_logic_corrSac_onset_045';
raster_data.train_data_logic_corrSac_vmax_045   = train_data_logic_corrSac_vmax_045';
raster_data.train_data_logic_corrSac_offset_045 = train_data_logic_corrSac_offset_045';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_090, 2)
raster_data.train_data_logic_SS_090 = train_data_logic_SS_090;
raster_data.train_data_logic_CS_090 = train_data_logic_CS_090;
raster_data.velocity_data_090       = velocity_data_090;
raster_data.train_data_logic_cue_present_090    = train_data_logic_cue_present_090;
raster_data.train_data_logic_primSac_onset_090  = train_data_logic_primSac_onset_090;
raster_data.train_data_logic_primSac_vmax_090   = train_data_logic_primSac_vmax_090;
raster_data.train_data_logic_primSac_offset_090 = train_data_logic_primSac_offset_090;
raster_data.train_data_logic_corrSac_onset_090  = train_data_logic_corrSac_onset_090;
raster_data.train_data_logic_corrSac_vmax_090   = train_data_logic_corrSac_vmax_090;
raster_data.train_data_logic_corrSac_offset_090 = train_data_logic_corrSac_offset_090;
else
raster_data.train_data_logic_SS_090 = train_data_logic_SS_090';
raster_data.train_data_logic_CS_090 = train_data_logic_CS_090';
raster_data.velocity_data_090       = velocity_data_090';
raster_data.train_data_logic_cue_present_090    = train_data_logic_cue_present_090';
raster_data.train_data_logic_primSac_onset_090  = train_data_logic_primSac_onset_090';
raster_data.train_data_logic_primSac_vmax_090   = train_data_logic_primSac_vmax_090';
raster_data.train_data_logic_primSac_offset_090 = train_data_logic_primSac_offset_090';
raster_data.train_data_logic_corrSac_onset_090  = train_data_logic_corrSac_onset_090';
raster_data.train_data_logic_corrSac_vmax_090   = train_data_logic_corrSac_vmax_090';
raster_data.train_data_logic_corrSac_offset_090 = train_data_logic_corrSac_offset_090';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_135, 2)
raster_data.train_data_logic_SS_135 = train_data_logic_SS_135;
raster_data.train_data_logic_CS_135 = train_data_logic_CS_135;
raster_data.velocity_data_135       = velocity_data_135;
raster_data.train_data_logic_cue_present_135    = train_data_logic_cue_present_135;
raster_data.train_data_logic_primSac_onset_135  = train_data_logic_primSac_onset_135;
raster_data.train_data_logic_primSac_vmax_135   = train_data_logic_primSac_vmax_135;
raster_data.train_data_logic_primSac_offset_135 = train_data_logic_primSac_offset_135;
raster_data.train_data_logic_corrSac_onset_135  = train_data_logic_corrSac_onset_135;
raster_data.train_data_logic_corrSac_vmax_135   = train_data_logic_corrSac_vmax_135;
raster_data.train_data_logic_corrSac_offset_135 = train_data_logic_corrSac_offset_135;
else
raster_data.train_data_logic_SS_135 = train_data_logic_SS_135';
raster_data.train_data_logic_CS_135 = train_data_logic_CS_135';
raster_data.velocity_data_135       = velocity_data_135';
raster_data.train_data_logic_cue_present_135    = train_data_logic_cue_present_135';
raster_data.train_data_logic_primSac_onset_135  = train_data_logic_primSac_onset_135';
raster_data.train_data_logic_primSac_vmax_135   = train_data_logic_primSac_vmax_135';
raster_data.train_data_logic_primSac_offset_135 = train_data_logic_primSac_offset_135';
raster_data.train_data_logic_corrSac_onset_135  = train_data_logic_corrSac_onset_135';
raster_data.train_data_logic_corrSac_vmax_135   = train_data_logic_corrSac_vmax_135';
raster_data.train_data_logic_corrSac_offset_135 = train_data_logic_corrSac_offset_135';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_180, 2)
raster_data.train_data_logic_SS_180 = train_data_logic_SS_180;
raster_data.train_data_logic_CS_180 = train_data_logic_CS_180;
raster_data.velocity_data_180       = velocity_data_180;
raster_data.train_data_logic_cue_present_180    = train_data_logic_cue_present_180;
raster_data.train_data_logic_primSac_onset_180  = train_data_logic_primSac_onset_180;
raster_data.train_data_logic_primSac_vmax_180   = train_data_logic_primSac_vmax_180;
raster_data.train_data_logic_primSac_offset_180 = train_data_logic_primSac_offset_180;
raster_data.train_data_logic_corrSac_onset_180  = train_data_logic_corrSac_onset_180;
raster_data.train_data_logic_corrSac_vmax_180   = train_data_logic_corrSac_vmax_180;
raster_data.train_data_logic_corrSac_offset_180 = train_data_logic_corrSac_offset_180;
else
raster_data.train_data_logic_SS_180 = train_data_logic_SS_180';
raster_data.train_data_logic_CS_180 = train_data_logic_CS_180';
raster_data.velocity_data_180       = velocity_data_180';
raster_data.train_data_logic_cue_present_180    = train_data_logic_cue_present_180';
raster_data.train_data_logic_primSac_onset_180  = train_data_logic_primSac_onset_180';
raster_data.train_data_logic_primSac_vmax_180   = train_data_logic_primSac_vmax_180';
raster_data.train_data_logic_primSac_offset_180 = train_data_logic_primSac_offset_180';
raster_data.train_data_logic_corrSac_onset_180  = train_data_logic_corrSac_onset_180';
raster_data.train_data_logic_corrSac_vmax_180   = train_data_logic_corrSac_vmax_180';
raster_data.train_data_logic_corrSac_offset_180 = train_data_logic_corrSac_offset_180';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_225, 2)
raster_data.train_data_logic_SS_225 = train_data_logic_SS_225;
raster_data.train_data_logic_CS_225 = train_data_logic_CS_225;
raster_data.velocity_data_225       = velocity_data_225;
raster_data.train_data_logic_cue_present_225    = train_data_logic_cue_present_225;
raster_data.train_data_logic_primSac_onset_225  = train_data_logic_primSac_onset_225;
raster_data.train_data_logic_primSac_vmax_225   = train_data_logic_primSac_vmax_225;
raster_data.train_data_logic_primSac_offset_225 = train_data_logic_primSac_offset_225;
raster_data.train_data_logic_corrSac_onset_225  = train_data_logic_corrSac_onset_225;
raster_data.train_data_logic_corrSac_vmax_225   = train_data_logic_corrSac_vmax_225;
raster_data.train_data_logic_corrSac_offset_225 = train_data_logic_corrSac_offset_225;
else
raster_data.train_data_logic_SS_225 = train_data_logic_SS_225';
raster_data.train_data_logic_CS_225 = train_data_logic_CS_225';
raster_data.velocity_data_225       = velocity_data_225';
raster_data.train_data_logic_cue_present_225    = train_data_logic_cue_present_225';
raster_data.train_data_logic_primSac_onset_225  = train_data_logic_primSac_onset_225';
raster_data.train_data_logic_primSac_vmax_225   = train_data_logic_primSac_vmax_225';
raster_data.train_data_logic_primSac_offset_225 = train_data_logic_primSac_offset_225';
raster_data.train_data_logic_corrSac_onset_225  = train_data_logic_corrSac_onset_225';
raster_data.train_data_logic_corrSac_vmax_225   = train_data_logic_corrSac_vmax_225';
raster_data.train_data_logic_corrSac_offset_225 = train_data_logic_corrSac_offset_225';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_270, 2)
raster_data.train_data_logic_SS_270 = train_data_logic_SS_270;
raster_data.train_data_logic_CS_270 = train_data_logic_CS_270;
raster_data.velocity_data_270       = velocity_data_270;
raster_data.train_data_logic_cue_present_270    = train_data_logic_cue_present_270;
raster_data.train_data_logic_primSac_onset_270  = train_data_logic_primSac_onset_270;
raster_data.train_data_logic_primSac_vmax_270   = train_data_logic_primSac_vmax_270;
raster_data.train_data_logic_primSac_offset_270 = train_data_logic_primSac_offset_270;
raster_data.train_data_logic_corrSac_onset_270  = train_data_logic_corrSac_onset_270;
raster_data.train_data_logic_corrSac_vmax_270   = train_data_logic_corrSac_vmax_270;
raster_data.train_data_logic_corrSac_offset_270 = train_data_logic_corrSac_offset_270;
else
raster_data.train_data_logic_SS_270 = train_data_logic_SS_270';
raster_data.train_data_logic_CS_270 = train_data_logic_CS_270';
raster_data.velocity_data_270       = velocity_data_270';
raster_data.train_data_logic_cue_present_270    = train_data_logic_cue_present_270';
raster_data.train_data_logic_primSac_onset_270  = train_data_logic_primSac_onset_270';
raster_data.train_data_logic_primSac_vmax_270   = train_data_logic_primSac_vmax_270';
raster_data.train_data_logic_primSac_offset_270 = train_data_logic_primSac_offset_270';
raster_data.train_data_logic_corrSac_onset_270  = train_data_logic_corrSac_onset_270';
raster_data.train_data_logic_corrSac_vmax_270   = train_data_logic_corrSac_vmax_270';
raster_data.train_data_logic_corrSac_offset_270 = train_data_logic_corrSac_offset_270';
end

if size(EPHYS_inds_event, 2) ==  size(train_data_logic_SS_315, 2)
raster_data.train_data_logic_SS_315 = train_data_logic_SS_315;
raster_data.train_data_logic_CS_315 = train_data_logic_CS_315;
raster_data.velocity_data_315       = velocity_data_315;
raster_data.train_data_logic_cue_present_315    = train_data_logic_cue_present_315;
raster_data.train_data_logic_primSac_onset_315  = train_data_logic_primSac_onset_315;
raster_data.train_data_logic_primSac_vmax_315   = train_data_logic_primSac_vmax_315;
raster_data.train_data_logic_primSac_offset_315 = train_data_logic_primSac_offset_315;
raster_data.train_data_logic_corrSac_onset_315  = train_data_logic_corrSac_onset_315;
raster_data.train_data_logic_corrSac_vmax_315   = train_data_logic_corrSac_vmax_315;
raster_data.train_data_logic_corrSac_offset_315 = train_data_logic_corrSac_offset_315;
else
raster_data.train_data_logic_SS_315 = train_data_logic_SS_315';
raster_data.train_data_logic_CS_315 = train_data_logic_CS_315';
raster_data.velocity_data_315       = velocity_data_315';
raster_data.train_data_logic_cue_present_315    = train_data_logic_cue_present_315';
raster_data.train_data_logic_primSac_onset_315  = train_data_logic_primSac_onset_315';
raster_data.train_data_logic_primSac_vmax_315   = train_data_logic_primSac_vmax_315';
raster_data.train_data_logic_primSac_offset_315 = train_data_logic_primSac_offset_315';
raster_data.train_data_logic_corrSac_onset_315  = train_data_logic_corrSac_onset_315';
raster_data.train_data_logic_corrSac_vmax_315   = train_data_logic_corrSac_vmax_315';
raster_data.train_data_logic_corrSac_offset_315 = train_data_logic_corrSac_offset_315';
end

end

%% function ESN_correlogram
function Corr_data = ESN_correlogram(SS_time, CS_time)
bin_size_time = 1e-3; % seconds
span_window_size = (1 / bin_size_time) * (100 / 1000);
span_window_size_half = round(span_window_size / 2);
inds_span = ((-span_window_size_half+1) : 1 : (span_window_size_half))';

if (~isempty(CS_time)) && (~isempty(SS_time))
    ch_time_min = min([SS_time(1) CS_time(1)]);
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max([SS_time(end) CS_time(end)]) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(CS_time))
    ch_time_min = min(  CS_time(1) );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max(  CS_time(end) ) + 2.0;
    
    CH__.CS_data.CS_time =  CS_time - ch_time_min;
    CH__.CS_data.CS_time = ESN_Round(CH__.CS_data.CS_time, bin_size_time);
    CH__.CS_data.CS_ind  = round(CH__.CS_data.CS_time .* (1/bin_size_time));
    CH__.CS_data.CS_ind( CH__.CS_data.CS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
elseif (~isempty(SS_time))
    ch_time_min = min( SS_time(1)  );
    ch_time_min = max([(ch_time_min-2.0) 0]);
    ch_time_max = max( SS_time(end)  ) + 2.0;
    
    CH__.SS_data.SS_time =  SS_time - ch_time_min;
    CH__.SS_data.SS_time = ESN_Round(CH__.SS_data.SS_time, bin_size_time);
    CH__.SS_data.SS_ind  = round(CH__.SS_data.SS_time .* (1/bin_size_time));
    CH__.SS_data.SS_ind( CH__.SS_data.SS_ind < 1 ) = 1;
    CH__.CH_data.ch_time_reconstruct = ch_time_min : bin_size_time : ch_time_max;
end

% SSxSS_AUTO
if (~isempty(SS_time))
    CH__.SS_data.SS_inds_reconstruct = repmat( CH__.SS_data.SS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.SS_data.SS_ind), 1);
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct < 1 ) = 1;
    CH__.SS_data.SS_inds_reconstruct( CH__.SS_data.SS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.SS_data.SS_inds_reconstruct );
    % SSxSS correlogram
    SSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    ss_inds_span     = repmat(inds_span(:)',     size(SS_time(:),1), 1);
    ss_bin_size_time = repmat(bin_size_time(:)', size(SS_time(:),1), 1);
else
    SSxSS_AUTO       = false(0, length(inds_span(:)'));
    ss_inds_span     = nan(0, length(inds_span(:)'));
    ss_bin_size_time = nan(0, 1);
end

% CSxSS_WITHIN
if (~isempty(CS_time)) && (~isempty(SS_time))
    CH__.CS_data.CS_inds_reconstruct = repmat( CH__.CS_data.CS_ind(:), 1, length(inds_span)) + repmat(inds_span(:)', length(CH__.CS_data.CS_ind), 1);
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct < 1 ) = 1;
    CH__.CS_data.CS_inds_reconstruct( CH__.CS_data.CS_inds_reconstruct > length( CH__.CH_data.ch_time_reconstruct ) ) = length( CH__.CH_data.ch_time_reconstruct );
    
    CH__.SS_data.SS_event_trace = false( size(CH__.CH_data.ch_time_reconstruct) );
    CH__.SS_data.SS_event_trace( CH__.SS_data.SS_ind ) = true ;
    CH__.SS_data.SS_event_trace( 1   ) = false;
    CH__.SS_data.SS_event_trace( end ) = false;
    
    CH__.SS_data.SS_event_reconstruct = CH__.SS_data.SS_event_trace( CH__.CS_data.CS_inds_reconstruct );
    % CSxSS correlogram
    CSxSS_AUTO       = CH__.SS_data.SS_event_reconstruct;
    cs_inds_span     = repmat(inds_span(:)',     size(CS_time(:),1), 1);
    cs_bin_size_time = repmat(bin_size_time(:)', size(CS_time(:),1), 1);
else
    CSxSS_AUTO       = false(0, length(inds_span(:)'));
    cs_inds_span     = nan(0, length(inds_span(:)'));
    cs_bin_size_time = nan(0, 1);
end

Corr_data = struct;
Corr_data.CS_inds_span     = cs_inds_span;
Corr_data.CS_bin_size_time = cs_bin_size_time;
Corr_data.SS_inds_span     = ss_inds_span;
Corr_data.SS_bin_size_time = ss_bin_size_time;
Corr_data.SS_SSxSS_AUTO    = SSxSS_AUTO;
Corr_data.CS_CSxSS_AUTO    = CSxSS_AUTO;
end