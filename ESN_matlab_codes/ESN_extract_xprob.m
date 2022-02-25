%% function extract_xprob
function [S1xS2_bool, output_span] = ESN_extract_xprob(spike1_bool, spike2_bool, sample_rate, bin_size, win_len_before, win_len_after)
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