%% function cambridge unit_extraction
function MAF_KS_unit_extraction(sess_path)

pat_SS = 'SS' + digitsPattern(1,2);
pat_CS = 'CS' + digitsPattern(1,2);
if ~strcmp(sess_path(end), filesep)
    sess_path = [sess_path filesep];
end

folders_ = strsplit(sess_path, filesep);
current_sess = folders_{end-1};
sess_name = regexprep(current_sess(3:end),'-','');
sess_meta_data = readtable([sess_path sess_name '.xls']);

rec_list = sess_meta_data.folder_name(logical...
    (sess_meta_data.ephys .* sess_meta_data.eye));

num_rec = length(rec_list);

for counter_rec = 1 : 1 : num_rec
    current_rec = rec_list{counter_rec};
    fprintf(['      ', num2str(counter_rec), ' / ' num2str(num_rec)...
        ' Analyzing rec ', current_rec '\n']);

    units_dir = [sess_path, current_rec, filesep, ...
            'analyzed_data', filesep, 'units' filesep];
    if exist(units_dir,'dir')
        rmdir(units_dir,'s');
    end

    path_to_sort = [sess_path, current_rec, filesep, ...
        'analyzed_data', filesep, 'sorted_data', filesep];
    path_to_data = [sess_path, current_rec, filesep, ...
        'raw_data', filesep];

    ch_file = dir([path_to_data,'*_CH1.continuous']);
    if isempty(ch_file)
        fprintf('\n No ephys!\n')
        continue;
    end
    [~, ch_time, ~] = load_open_ephys_data(...
        [path_to_data,ch_file(1).name]);
    t_start = ch_time(1);
    n_samples = length(ch_time);
    sp = loadKSdir(path_to_sort);
    sp_inds = sp.sp_info.st;
    sp_clust = sp.sp_info.clu;
    sample_rate = sp.sample_rate;
    cl_waveform = sp.waveform;

    ch_num = sp.clu_info.ch;
    id = sp.clu_info.id;
    [~, ind_good] = ismember(sp.clu_info.group, 'good', 'rows');
    good_ids = id(ind_good == 1);
    neurontype = sp.clu_info.neurontype;
    rec_name = rec_list{counter_rec};

    for counter_gid = 1:numel(good_ids)
        current_id = good_ids(counter_gid);
        label = deblank(neurontype(id == current_id,:));
        if startsWith(label,pat_SS)
            ss_id = current_id;
            pc_num = str2double(erase(label,'SS'));
            [~, ind_cs] = ismember(neurontype,...
                ['CS' num2str(pc_num)],'rows');
            cs_id = id(ind_cs == 1);
            ch_id = cs_id;
            type = 'PC';
        elseif startsWith(label,pat_CS)
            continue;
        elseif strcmp(label,'CS')
            ss_id = nan;
            cs_id = current_id;
            ch_id = cs_id;
            type = 'CS';
        else
            ss_id = current_id;
            cs_id = nan;
            ch_id = ss_id;
            type = label;
        end

        ss_index_int = sp_inds(sp_clust == ss_id);
        cs_index_int = sp_inds(sp_clust == cs_id);
        ch_ = ch_num(id == ch_id);

        time_ind = 1:n_samples;
        ss_index = ismember(time_ind,ss_index_int)';
        cs_index = ismember(time_ind,cs_index_int)';

        ss_index = resolve_auto_conflicts(ss_index,sample_rate,.5e-3);

        if sum(cs_index) > 0
            cs_index = resolve_auto_conflicts(cs_index,sample_rate,5e-3);
            ss_index = resolve_cs_ss_conflicts(ss_index,cs_index,sample_rate,.5e-3,.5e-3);
        end

        waveform.ch_map.x = sp.xcoords;
        waveform.ch_map.y = sp.ycoords;
        waveform.ch_map.map = sp.ch_map;

        if sum(ss_index) > 0
            ss_wave = squeeze(cl_waveform.mean_wf(id == ss_id,:,:))';
            ss_wave_ch = squeeze(cl_waveform.ch_num(id == ss_id,:))';
            ss_peak = sp.sp_info.amp(sp_clust == ss_id);
        end
        
        if sum(cs_index) > 0
            cs_wave = squeeze(cl_waveform.mean_wf(id == cs_id,:,:))';
            cs_wave_ch = squeeze(cl_waveform.ch_num(id == cs_id,:))';
            cs_peak = sp.sp_info.amp(sp_clust == cs_id);
        end

        if sum(ss_index) > 0
            waveform.ss_wave    = ss_wave;
            waveform.ss_wave_ch = ss_wave_ch;
            waveform.ss_peak    = ss_peak;
        else
            waveform.ss_wave    = nan(size(cs_wave));
            waveform.ss_wave_ch = cs_wave_ch;
            waveform.ss_peak    = [];
        end

        if sum(cs_index) > 0
            waveform.cs_wave    = cs_wave;
            waveform.cs_wave_ch = cs_wave_ch;
            waveform.cs_peak    = cs_peak;
        else
            waveform.cs_wave    = nan(size(ss_wave));
            waveform.cs_wave_ch = ss_wave_ch;
            waveform.cs_peak    = [];
        end

        sess_name = regexprep(rec_name(3:end),'-','');
        unit_name = [sess_name '_' num2str(ch_+1,'%.2i') '_' ...
            num2str(current_id,'%.3i')];
        out_dir = [units_dir unit_name filesep];
        mkdir(out_dir);
        save([out_dir unit_name '_sorted_KSA'],'waveform','t_start',...
            'sample_rate','ss_index','cs_index','type','-v7.3');
        close all;
        MAF_plot_cell_summary([out_dir unit_name '_sorted_KSA']);

    end

end
end

%% function solve SS-SS Conflicts
function ss_index_new = resolve_auto_conflicts(ss_index,sample_rate,win_look_around)
window_len = floor(win_look_around * sample_rate);
ss_index_ = ss_index;
ss_index_int_ = find(ss_index_);
for counter_ss = 1:length(ss_index_int_)
    ss_index_local_ = ss_index_int_(counter_ss);
    if ss_index_local_ < window_len
        ss_index_(ss_index_local_) = 0;
        continue
    end
    if ss_index_local_ > (length(ss_index_) - window_len)
        ss_index_(ss_index_local_) = 0;
        continue
    end
    search_win_inds = ss_index_local_-window_len:ss_index_local_+window_len;
    ss_search_win_bool = ss_index_(search_win_inds);
    ss_search_win_int  = find(ss_search_win_bool);
    if length(ss_search_win_int) < 2
        continue
    end
    if length(ss_search_win_int) > 1
        valid_ind = ss_search_win_int(1);
        ss_search_win_bool = zeros(size(search_win_inds));
        ss_search_win_bool(valid_ind) = 1;
        ss_index_(search_win_inds) = ss_search_win_bool;
    end
end
ss_index_new = ss_index_;
end

%% function solve SS-CS Conflicts
function ss_index_new = resolve_cs_ss_conflicts(ss_index,cs_index,sample_rate, win_look_before, win_look_after)
window_len_back = floor(win_look_before * sample_rate);
window_len_front = floor(win_look_after * sample_rate);
cs_index_int = find(cs_index);
ss_index_ = ss_index;
for counter_cs = 1:length(cs_index_int)
    cs_index_local = cs_index_int(counter_cs);
    search_win_inds = cs_index_local-window_len_back:cs_index_local+window_len_front;
    ss_search_win_bool = ss_index_(search_win_inds);
    ss_search_win_int  = find(ss_search_win_bool);
    if ~isempty(ss_search_win_int)
        ss_ind_invalid = ss_search_win_int + cs_index_local - window_len_back;
        ss_index_(ss_ind_invalid-1) = 0;
    end
end
ss_index_new = ss_index_;
end
