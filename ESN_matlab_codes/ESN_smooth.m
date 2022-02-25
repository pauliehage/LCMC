function smooth_data_ = ESN_smooth(data_, dim)
% smooth data using 2nd order Savitzky-Golay filter with 21 points
% if data_ is a matrix, the method will smooth each column by default or smooth along dim.
if nargin < 2
    dim = 1;
end
if size(data_, 1) == 1
    smooth_data_ = reshape(smooth_data(data_), 1, []);
elseif size(data_, 2) == 1
    smooth_data_ = reshape(smooth_data(data_), [], 1);
else
    smooth_data_ = nan(size(data_));
    if dim == 1
        % smooth columns
        for counter_col = 1 : size(data_, 2)
            smooth_data_(:, counter_col) = reshape(smooth_data(data_(:, counter_col)), [], 1);
        end
    elseif dim == 2
        % smooth rows
        for counter_row = 1 : size(data_, 1)
            smooth_data_(counter_row, :) = reshape(smooth_data(data_(counter_row, :)), 1, []);
        end
    end
    
end
end

function smooth_data_ = smooth_data(data_)
% method = 'moving';  % Moving average. A lowpass filter with filter coefficients equal to the reciprocal of the span.
% method = 'lowess';  % Local regression using weighted linear least squares and a 1st degree polynomial model.
% method = 'loess';   % Local regression using weighted linear least squares and a 2nd degree polynomial model.
% method = 'sgolay';  % Savitzky-Golay filter. A generalized moving average with filter coefficients determined by an unweighted linear least-squares regression and a polynomial model of specified degree (default is 2). The method can accept nonuniform predictor data.
% method = 'rlowess'; % A robust version of 'lowess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% method = 'rloess';  % A robust version of 'loess' that assigns lower weight to outliers in the regression. The method assigns zero weight to data outside six mean absolute deviations.
% smooth_data_ = smooth(data_, method);
%
order = 2;
num_points = 31;
smooth_data_ = smooth(data_, num_points, 'sgolay', order);
%}

% filter params
%{
if sum(isnan(data_))
    smooth_data_ = data_;
    return;
end
order = 3;
sampling_freq = 1000.0;
cutoff_freq = 50.0;
[b_butter,a_butter] = butter(order,(cutoff_freq/(sampling_freq/2)), 'low');
smooth_data_ = filtfilt(b_butter,a_butter,data_);
%}
end