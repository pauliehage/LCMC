function DB = MAF_extract_cell_data(DB)

cs_index = DB.cs_index;
ss_index = DB.ss_index;
sample_rate = DB.sample_rate;

GLOBAL_XPROB_SS_BINSIZE    = 1e-3;
GLOBAL_XPROB_SS_BEFORE     = 5e-2;
GLOBAL_XPROB_SS_AFTER      = 5e-2;
GLOBAL_XPROB_CS_BINSIZE    = 1e-3;
GLOBAL_XPROB_CS_BEFORE     = 5e-2;
GLOBAL_XPROB_CS_AFTER      = 5e-2;

ss_isi = extract_isi(ss_index, sample_rate);
cs_isi = extract_isi(cs_index, sample_rate);
[ss_xprob, ss_xprob_span] = extract_xprob(ss_index, ss_index, sample_rate, ...
    GLOBAL_XPROB_SS_BINSIZE, ...
    GLOBAL_XPROB_SS_BEFORE, ...
    GLOBAL_XPROB_SS_AFTER);
win_len_before_int = round(double(GLOBAL_XPROB_SS_BEFORE) ...
    / double(GLOBAL_XPROB_SS_BINSIZE));
ss_xprob(:,win_len_before_int+1) = NaN;
[cs_xprob, cs_xprob_span] = extract_xprob(cs_index, ss_index, sample_rate, ...
    GLOBAL_XPROB_CS_BINSIZE, ...
    GLOBAL_XPROB_CS_BEFORE, ...
    GLOBAL_XPROB_CS_AFTER);

% store data
DB.ss_isi = ss_isi;
DB.cs_isi = cs_isi;
DB.ss_xprob = ss_xprob;
DB.ss_xprob_span = ss_xprob_span;
DB.cs_xprob = cs_xprob;
DB.cs_xprob_span = cs_xprob_span;
end

%% function extract_isi
function inter_spike_interval = extract_isi(index_bool, sample_rate)
index_bool = logical(index_bool); index_bool(1) = false; index_bool(end) = false;
if sum(index_bool) < 2
    inter_spike_interval = [];
    return;
end
index_value = find(index_bool);
inter_spike_interval = diff(index_value) / double(sample_rate);
inter_spike_interval = [inter_spike_interval(:); inter_spike_interval(end)];
end

%% function extract_xprob
function [S1xS2_bool, output_span] = extract_xprob(spike1_bool, spike2_bool, sample_rate, bin_size, win_len_before, win_len_after)
spike1_bool = logical(spike1_bool); spike1_bool(1) = false; spike1_bool(end) = false;
spike2_bool = logical(spike2_bool); spike2_bool(1) = false; spike2_bool(end) = false;
spike1_time = find(spike1_bool) / double(sample_rate);
spike1_int = round(spike1_time/double(bin_size));
spike2_time = find(spike2_bool) / double(sample_rate);
spike2_index = round(spike2_time/double(bin_size));
spike2_bool_size = round(double(length(spike1_bool)) / double(sample_rate) / double(bin_size));
spike2_bool = zeros(spike2_bool_size,1);
spike2_index(spike2_index<1) = 1;
spike2_index(spike2_index>spike2_bool_size) = spike2_bool_size;
spike2_bool(spike2_index) = 1;
win_len_before_int = round(double(win_len_before) / double(bin_size));
win_len_after_int = round(double(win_len_after) / double(bin_size));
span_int = -win_len_before_int : 1 : win_len_after_int;
num_row = length(spike1_int);
num_col = length(span_int);
spike1_int = repmat(spike1_int(:), 1, num_col);
span_int = repmat(span_int(:)', num_row, 1);
ind = spike1_int + span_int;
ind(ind<1) = 1;
ind(ind>length(spike2_bool)) = length(spike2_bool);
S1xS2_bool = spike2_bool(ind);
output_span = span_int * double(bin_size);
if sum(spike1_bool)==1
    S1xS2_bool  = reshape(S1xS2_bool, 1, []);
    output_span = reshape(output_span, 1, []);
end
end