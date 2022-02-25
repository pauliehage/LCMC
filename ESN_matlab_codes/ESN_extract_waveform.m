%% function extract_waveform
function [waveform, span] = ESN_extract_waveform(data, spike_bool, sample_rate, win_len_before, win_len_after)
spike_bool = logical(spike_bool); spike_bool(1) = false; spike_bool(end) = false;
if sum(spike_bool) < 1
    waveform = [];
    span = [];
    return;
end
spike_int = find(spike_bool);
win_len_before_int = round(double(win_len_before) * double(sample_rate));
win_len_after_int  = round(double(win_len_after)  * double(sample_rate));
span_int = (-win_len_before_int : 1 : win_len_after_int);
num_row = length(spike_int);
num_col = length(span_int);
spike_int = repmat(spike_int(:), 1, num_col);
span_int  = repmat(span_int(:)', num_row, 1);

ind = spike_int + span_int;
ind(ind<1) = 1;
ind(ind>length(data)) = length(data);
waveform = data(ind);
span = span_int / double(sample_rate);
if sum(spike_bool)==1
    waveform = reshape(waveform, 1, []);
    span     = reshape(span, 1, []);
end
end
