function [] = PGH_write_repeated_measures_csv(data, filename, varargin)
% WRITE_REPEATED_MEASURES_CSV
%   [] = WRITE_REPEATED_MEASURES_CSV(data, filename, varargin) writes
%   the contents of the input structure data to a csv file whose filename
%   is specified by the filename parameter.
%   
%   In the simplest case, data is an NxM data array containing M repeated
%   measures for N subjects. The id's for each subject are ordered from 1
%   to M as the first column in the output file.
%
%   If `data' is a cell array, we write each of the data arrays inside the
%   cell array as a different group (group numbering starts from 1 to the
%   length of the cell array).
%
%   Additional key-value pairs can be specified:
%   -'long' Use a long format. Each of the repeated measurements is
%   specified as an additional row in the output file. By default, we store
%   the output in short form. That is, each independent sample (i.e.
%   subject) is stored as a single row. 
%   -'measurement_axis': Specify the name of the repeated measures. That
%   is, in short mode, each of the repeated measures is defined by a name
%   which defaults to 'measure_#' (where # is the number of the repeated
%   measure). This parameter changes the name of that measurement label
%   -'id' Specifies the ID for each of the measurements, N. If data is a
%   cell array, this is also expected to be a cell array of ids.
%

args = getopt(varargin, {'long', {'measurement_axis=' []}, {'id=', []}});

% Create the stats directory
temp_path = mfilename('fullpath');
name_ind = strfind(temp_path,mfilename);
path = temp_path(1:name_ind-1);
stats_directory = strcat(path,'stats');
if ~exist(stats_directory, 'dir')
    mkdir(stats_directory);
end

filename = strcat(stats_directory,'\', filename);

% Make sure the ids are specified properly
if isempty(args.id)
    if ~iscell(data)
        args.id = 1:size(data, 1);
    else
        args.id = cellfun(@(x) 1:size(x, 1), data, 'UniformOutput', 0);
        for i = 1:length(args.id)
            args.id{i} = args.id{i} + size(data{i}, 1) * (i - 1);
        end
    end
else
    if iscell(data) && ~iscell(args.id)
        error('Invalid ids specified. Expected a cell array');
    end
end

% Make sure repeated dimensions are the same
if iscell(data)
    dims = size(data{1});
    for i = 1:numel(data)
        if all(size(data{i}) ~= dims)
            error('Size of each cell of data are not consistent')
        end
    end
    cellDims = ndims(data);
end

if ~isempty(args.measurement_axis) && ((~iscell(data) && length(args.measurement_axis) ~= size(data, 2)) || iscell(data) && length(args.measurement_axis) ~= size(data{1}, 2))
    error('Invalid size of measurement axis')
end

% Open the file for writing
fp = fopen(filename, 'w');

fprintf(fp, '"subject"');

if (iscell(args.id) && isnumeric(args.id{1})) || (~iscell(args.id) && isnumeric(args.id))
    format_string = '%d';
else
    format_string = '"%s"';
end

if iscell(data)
    if cellDims > 1
        fprintf(fp, ' "phase" "errDir"');
        format_string = strcat(format_string, ' %d %d');        
    else
        fprintf(fp, ' "group"');
        format_string = strcat(format_string, ' %d');
    end
end

if args.long
    fprintf(fp, ' "trial" "outcome"');
    format_string = strcat(format_string, ' %d %f');
else % Use short (row) format
    if iscell(data)
        temp = data{1};
    else
        temp = data;
    end
    for i = 1:size(temp, 2)
        if isempty(args.measurement_axis)
            fprintf(fp, [' "measure_' num2str(i) '"']);
        else
            fprintf(fp, [' "measure_' num2str(args.measurement_axis(i)) '"']);
        end
        format_string = strcat(format_string, ' %f');
    end
end
fprintf(fp, '\r\n');
format_string = strcat(format_string, '\r\n');


if args.long
    if iscell(data)
        for i = 1:size(data,1)
            for ii = 1:size(data,2)
                temp = data{i,ii};
                for j = 1:size(temp, 1)
                    for k = 1:size(temp, 2)
                        fprintf(fp, format_string, args.id{i,ii}(j), i, ii, k, temp(j, k));
                    end
                end
            end
        end
    else
        for i = 1:size(data, 1)
            for j = 1:size(data, 2)
                fprintf(fp, format_string, args.id(i), j, data(i, j));
            end
        end
    end
else
    if iscell(data)
        for i = 1:length(data)
            temp = data{i};
            for j = 1:size(temp, 1)
                if isnumeric(args.id{i})
                    fprintf(fp, format_string, args.id{i}(j), i, temp(j, :));
                else
                    fprintf(fp, format_string, args.id{i}{j}, i, temp(j, :));
                end
            end
        end
    else
        % Write out the data
        for i = 1:size(data, 1)
            if isnumeric(args.id)
                fprintf(fp, format_string, args.id(i), data(i, :));
            else
                fprintf(fp, format_string, args.id{i}, data(i, :));
            end
        end
    end
end

fclose(fp);
