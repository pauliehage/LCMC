%% function psort_extract_unit
function MAF_psort_unit_extraction(sess_path)

if ~strcmp(sess_path(end), filesep)
    sess_path = [sess_path filesep];
end

wave_span = [2e-3,4e-3];

folders_ = strsplit(sess_path,filesep);
current_sess = folders_{end-1};
sess_name = regexprep(current_sess(3:end),'-','');
sess_meta_data = readtable([sess_path sess_name '.xls']);

rec_list = sess_meta_data.folder_name(logical...
    (sess_meta_data.ephys .* sess_meta_data.eye));

elec_list = sess_meta_data.elec(logical...
    (sess_meta_data.ephys .* sess_meta_data.eye));
elec = elec_list(1);

num_rec = length(rec_list);

for counter_rec = 1 : 1 : num_rec
    current_rec = rec_list{counter_rec};
    fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec)...
        ' Analyzing rec ', current_rec ' \n']);

    path_to_analyzed = [sess_path, current_rec, filesep, ...
        'analyzed_data', filesep];
    path_to_sort = [path_to_analyzed, 'sorted_data', filesep];
    path_to_data = [sess_path, current_rec, filesep, ...
        'raw_data', filesep];

    ch_file = dir([path_to_data,'*_CH*.continuous']);
    n_ch = numel(ch_file);
    if n_ch == 0
        fprintf('\n No ephys!\n')
        continue;
    end

    % specify geom
    if elec == 3
        x = 3*[sqrt(3),sqrt(3),sqrt(3)/2,3*sqrt(3)/2]';
        y = 6*[3,2,3/2,3/2]';
        ch_map = [1,2,3,4]-1;
        if not(n_ch == 4)
            error('inconsistent channel numbers!');
        end
    elseif elec == 4
        x = 3*[0, sqrt(3), sqrt(3)/2, sqrt(3), 3*sqrt(3)/2, sqrt(3), 2*sqrt(3)]';
        y = 6*[3, 0, 3/2, 2, 3/2, 3, 3]';
        ch_map = [1, 2, 3, 4, 5, 6, 7]-1;
        if not(n_ch == 7)
            error('inconsistent channel numbers!');
        end
    end

    units_list = dir([path_to_sort '*.psort']);
    units_list = {units_list.name};

    flag_ch = zeros(n_ch,1);
    units_struct = struct;
    num_units = numel(units_list);

    ch_data = [];

    for counter_unit = 1:num_units
        current_unit = units_list{counter_unit};
        DB = Psort_read_psort([path_to_sort current_unit]);
        ss_index = logical(DB.topLevel_data.ss_index);
        cs_index = logical(DB.topLevel_data.cs_index);

        if sum(cs_index) > 10 && sum(ss_index) > 0
            type = 'PC';
        elseif sum(cs_index) > 10
            type = 'CS';
        else
            cs_index(cs_index == 1) = 0;
            type = '';
        end

        ch = str2double(current_unit(15:16));

        if counter_unit == 1
            ch_data = zeros(n_ch,length(ss_index));
            ch_time = DB.topLevel_data.ch_time;
            t_start = ch_time(1);
            sample_rate = double(DB.topLevel_data.sample_rate);
        end

        ch_data(ch,:) = DB.topLevel_data.ch_data;
        flag_ch(ch) = 1;

        units_struct(counter_unit).ss_index = ss_index;
        units_struct(counter_unit).cs_index = cs_index;
        units_struct(counter_unit).type = type;

    end

    median_file = dir([path_to_analyzed '*_median.h5']);
    median = h5read([path_to_analyzed median_file.name],'/ch_data');

    for counter_ch = 1:n_ch
        if flag(counter_ch) == 1
            continue;
        end
        ch_file = dir([path_to_data ...
            '*_CH' num2str(counter_ch*4-3) '.continuous']);
        if isempty(ch_file)
            ch_file = dir([path_to_data ...
                '*_CH' num2str(counter_ch*4) '.continuous']);
        end
        [ch_data_, ~, ~] = load_open_ephys_data(...
            [path_to_data,ch_file(1).name]);
        ch_data(counter_ch,:) = ch_data_ - median;
    end

    for counter_unit = 1:num_units

        current_unit = units_list{counter_unit};

        fprintf(['                      ', num2str(counter_unit),...
            ' / ' num2str(num_units) ' Analyzing unit ',...
            current_unit ' \n']);

        ss_index = units_struct(counter_unit).ss_index;
        cs_index = units_struct(counter_unit).cs_index;
        type = units_struct(counter_unit).type;

        waveform.ch_map.x = x;
        waveform.ch_map.y = y;
        waveform.ch_map.map = ch_map;

        waveform_length = round(double(wave_span(1)) *...
            double(sample_rate)) + ...
            round(double(wave_span(2)) * double(sample_rate)) + 1;

        ch_ = str2double(current_unit(15:16));

        waveform.ss_wave_ch = ch_map;
        waveform.ss_wave = nan(n_ch,waveform_length);
        waveform.ss_peak = ch_data(ch_,logical(ss_index));

        waveform.cs_wave_ch = ch_map;
        waveform.cs_wave = nan(n_ch,waveform_length);
        waveform.cs_peak = ch_data(ch_,logical(cs_index));

        if sum(ss_index) > 0
            for counter_ch = 1:n_ch
                [waveform_, ~] = ESN_extract_waveform(...
                    ch_data(counter_ch,:), ss_index, sample_rate,...
                    wave_span(1), wave_span(2));
                waveform.ss_wave(counter_ch,:) = mean(waveform_);
            end
        end

        if sum(cs_index) > 0
            for counter_ch = 1:n_ch
                [waveform_, ~] = ESN_extract_waveform(...
                    ch_data(counter_ch,:), cs_index, sample_rate,...
                    wave_span(1), wave_span(2));
                waveform.cs_wave(counter_ch,:) = mean(waveform_);
            end
        end

        user = current_unit(end-8:end-6);

        current_id = counter_unit;
        sess_name = regexprep(current_rec(3:end),'-','');
        unit_name = [sess_name '_' num2str(ch_,'%.2i') '_' num2str(current_id,'%.3i')];
        out_dir = [path_to_analyzed,...
            'units', filesep unit_name filesep];
        mkdir(out_dir);
        save([out_dir unit_name '_sorted_' user],'waveform','t_start',...
            'sample_rate','ss_index','cs_index','type','-v7.3');
        close all;
        MAF_plot_cell_summary([out_dir unit_name '_sorted_' user])

    end

end

plot_rec_sort_summary(sess_path)
end

%% function extract session summary plot
function plot_rec_sort_summary(sess_path)

GLOBAL_XPROB_SS_BINSIZE    = 1e-3;
GLOBAL_XPROB_SS_BEFORE     = 5e-2;
GLOBAL_XPROB_SS_AFTER      = 5e-2;
GLOBAL_XPROB_CS_BINSIZE    = 1e-3;
GLOBAL_XPROB_CS_BEFORE     = 5e-2;
GLOBAL_XPROB_CS_AFTER      = 5e-2;

if ~strcmp(sess_path(end), filesep)
    sess_path = [sess_path filesep];
end

folders_ = strsplit(sess_path,filesep);
current_sess = folders_{end-1};
sess_name = regexprep(current_sess(3:end),'-','');

sess_meta_data = readtable([sess_path sess_name '.xls']);
rec_list = sess_meta_data.folder_name(logical...
    (sess_meta_data.ephys .* sess_meta_data.eye));

num_rec = length(rec_list);

for counter_rec = 1 : 1 : num_rec
    current_rec = rec_list{counter_rec};
    fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec) ' Analyzing rec ', current_rec ' \n']);

    units_path = [sess_path, current_rec, filesep, ...
        'analyzed_data', filesep 'units' filesep];

    sess_name = regexprep(current_rec(3:end),'-','');

    units_list = dir([units_path , sess_name '*']);
    units_list([units_list.isdir] == 0) = [];
    units_list = {units_list.name};

    num_units = length(units_list);
    units = struct;

    close all;
    fig_handle = figure('Units','normalized','Position',[0,0,1,1]);
    clf(fig_handle);
    num_row = num_units;
    num_col = num_units+1;
    sgtitle(fig_handle, [sess_name, ' (i,j): P(j|i=0)'],...
        'Interpreter','none');

    for counter_unit = 1:num_units
        current_unit = units_list{counter_unit};
        sorted_file = dir([units_path current_unit filesep ...
            '*_sorted_*.mat']);
        load([units_path current_unit filesep ...
            sorted_file(1).name],'cs_index','ss_index','t_start',...
            'type','waveform','sample_rate');
        units(counter_unit).name = sorted_file(1).name(15:end-4);
        units(counter_unit).ss_index = ss_index;
        units(counter_unit).cs_index = cs_index;
        units(counter_unit).t_start = t_start;
        units(counter_unit).type = type;
        units(counter_unit).waveform = waveform;
        units(counter_unit).sample_rate = sample_rate;
        subplot(num_row, num_col, (counter_unit-1)*num_col+1);
        plot_waveform(waveform,sample_rate)
        title(units(counter_unit).name,'Interpreter','none');
    end

    for counter_row = 1:num_units
        subplot(num_row,num_col,(counter_row-1)*num_col+counter_row+1)
        if counter_row == 1
            title(units(counter_row).name,'Interpreter','none');
        end
        hold on;
        unit_row = units(counter_row);
        ss_index = unit_row.ss_index;
        cs_index = unit_row.cs_index;
        flag_ss_row = sum(ss_index)>0;
        flag_cs_row = sum(cs_index)>0;
        if flag_ss_row
            [ss_xprob, ss_xprob_span] = ESN_extract_xprob(ss_index,...
                ss_index, sample_rate, ...
                GLOBAL_XPROB_SS_BINSIZE, ...
                GLOBAL_XPROB_SS_BEFORE, ...
                GLOBAL_XPROB_SS_AFTER);
            win_len_before_int = round(double(GLOBAL_XPROB_SS_BEFORE) ...
                / double(GLOBAL_XPROB_SS_BINSIZE));
            ss_xprob(:,win_len_before_int+1) = NaN;
            % plot
            ss_xprob_span_mean = mean(ss_xprob_span);
            ss_xprob_mean = mean(ss_xprob);
            plot(ss_xprob_span_mean*1000, ss_xprob_mean, '-b',...
                'linewidth', 2);
        elseif flag_cs_row
            [cs_xprob, cs_xprob_span] = ESN_extract_xprob(cs_index, cs_index, sample_rate, ...
                GLOBAL_XPROB_SS_BINSIZE, ...
                GLOBAL_XPROB_SS_BEFORE, ...
                GLOBAL_XPROB_SS_AFTER);
            win_len_before_int = round(double(GLOBAL_XPROB_SS_BEFORE) ...
                / double(GLOBAL_XPROB_SS_BINSIZE));
            cs_xprob(:,win_len_before_int+1) = NaN;
            % plot
            cs_xprob_span_mean = mean(cs_xprob_span);
            cs_xprob_mean = mean(cs_xprob);
            plot(cs_xprob_span_mean*1000, cs_xprob_mean, '-r',...
                'linewidth', 2);
        end
        if flag_cs_row && flag_ss_row
            [cs_xprob, cs_xprob_span] = ESN_extract_xprob(cs_index, ss_index, sample_rate, ...
                GLOBAL_XPROB_CS_BINSIZE, ...
                GLOBAL_XPROB_CS_BEFORE, ...
                GLOBAL_XPROB_CS_AFTER);
            cs_xprob_span_mean = mean(cs_xprob_span);
            cs_xprob_mean = mean(cs_xprob);
            plot(cs_xprob_span_mean*1000, cs_xprob_mean, '-r',...
                'linewidth', 2);
        end

        for counter_col = 1:num_units
            if counter_col == counter_row
                continue;
            end
            subplot(num_row,num_col,(counter_row-1)*num_col+counter_col+1)
            if counter_row == 1
                title(units(counter_col).name,'Interpreter','none');
            end
            hold on;
            unit_col = units(counter_col);
            flag_ss_col = sum(unit_col.ss_index)>0;
            flag_cs_col = sum(unit_col.cs_index)>0;
            if flag_ss_row
                spike_row = unit_row.ss_index;
                if flag_ss_col
                    spike_col = unit_col.ss_index;
                    [xprob, xprob_span] = ESN_extract_xprob(spike_row, spike_col, sample_rate, ...
                        GLOBAL_XPROB_CS_BINSIZE, ...
                        GLOBAL_XPROB_CS_BEFORE, ...
                        GLOBAL_XPROB_CS_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    plot(xprob_span_mean*1000, xprob_mean, '-b',...
                        'linewidth', 2);
                end
                if flag_cs_col
                    spike_col = unit_col.cs_index;
                    [xprob, xprob_span] = ESN_extract_xprob(spike_row,...
                        spike_col, sample_rate, ...
                        GLOBAL_XPROB_CS_BINSIZE, ...
                        GLOBAL_XPROB_CS_BEFORE, ...
                        GLOBAL_XPROB_CS_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    plot(xprob_span_mean*1000, xprob_mean, '-r',...
                        'linewidth', 2);
                end
            end
            if flag_cs_row
                spike_row = unit_row.cs_index;
                if flag_ss_col
                    spike_col = unit_col.ss_index;
                    [xprob, xprob_span] = ESN_extract_xprob(spike_row,...
                        spike_col, sample_rate, ...
                        GLOBAL_XPROB_CS_BINSIZE, ...
                        GLOBAL_XPROB_CS_BEFORE, ...
                        GLOBAL_XPROB_CS_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    plot(xprob_span_mean*1000, xprob_mean, '-m',...
                        'linewidth', 2);
                end
                if flag_cs_col
                    spike_col = unit_col.cs_index;
                    [xprob, xprob_span] = ESN_extract_xprob(spike_row,...
                        spike_col, sample_rate, ...
                        GLOBAL_XPROB_CS_BINSIZE, ...
                        GLOBAL_XPROB_CS_BEFORE, ...
                        GLOBAL_XPROB_CS_AFTER);
                    xprob_span_mean = mean(xprob_span);
                    xprob_mean = mean(xprob);
                    plot(xprob_span_mean*1000, xprob_mean, '-k',...
                        'linewidth', 2);
                end
            end
        end
    end
    ESN_Beautify_Plot(fig_handle, [15,9])
    fig_path = [sess_path, current_rec, filesep, ...
        'analyzed_figs' filesep];
    mkdir(fig_path);
    saveas(fig_handle,[fig_path sess_name '_sort_summary.pdf'],'pdf');
end
end

%% function plot_waveform
function plot_waveform(waveform,sample_rate)
hold on
ss_wave = waveform.ss_wave;
cs_wave = waveform.cs_wave;

x_ = waveform.ch_map.x * 4;
y_ = waveform.ch_map.y * 100;

if not(isnan(cs_wave))
    n_ch = size(cs_wave,1);
    n_sig = length(cs_wave);
    wave_ch = waveform.cs_wave_ch;
elseif not(isnan(ss_wave))
    n_ch = size(ss_wave,1);
    n_sig = length(ss_wave);
    wave_ch = waveform.ss_wave_ch;
end

x = x_(wave_ch+1);
y = y_(wave_ch+1);

span_ind = (0:n_sig-1)/sample_rate;
span_group_ = repmat([span_ind,nan],n_ch,1);
span_group = reshape((span_group_*1e3+x)',1,n_ch*(n_sig+1));

if not(isnan(ss_wave))
    ss_wave_ = [ss_wave, nan(n_ch,1)];
    ss_wave_group = reshape((ss_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, ss_wave_group, '-b', 'linewidth', 2)
end

if not(isnan(cs_wave))
    cs_wave_ = [cs_wave, nan(n_ch,1)];
    cs_wave_group = reshape((cs_wave_+y)',1,n_ch*(n_sig+1));
    plot(span_group, cs_wave_group, '-r', 'linewidth', 2);
end

ch_map = arrayfun(@num2str, wave_ch+1,'UniformOutput', 0);

text(x-1,y,ch_map)

axis off

plot([0,1],[0,0]+min(y)-100,'black','LineWidth',1);

end