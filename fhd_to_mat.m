function fhd_to_mat(folder_path)
% get the folder path and convert all the fhd files in that folder to mat files

% if there is no inputs, then set folder_path to pwd
if nargin < 1
    folder_path = pwd;
end
% add '/' to the end of folder_path
if ~strcmp(folder_path(end), '/')
    folder_path = [folder_path '/'];
end
% get the list of fhd files to analyze
FILES_FHD = dir([folder_path '*.fhd']);

for counter_h5Files = 1 : length(FILES_FHD)
    file_fhd = FILES_FHD(counter_h5Files);
    [~,file_name,~] = fileparts(file_fhd.name);
    file_path_name = [folder_path file_fhd.name];
    
    fprintf([num2str(counter_h5Files) '. Loading: ' file_name ' ... '])
    
    clearvars('data');
    data = fhd_load(file_path_name);
    save([folder_path '/' file_name '.mat'], 'data' );
    
    fprintf(' --> Completed.\n')
end
fprintf('ALL Completed. \n')
end

function [data] = fhd_load(filename, varargin)
% Load an FHD (Fast Hierarchical Data) file into a structure.
%
%   [data] = fhd_load(filename) loads all data from an FHD data file.
%

if exist('getopt', 'file')
    args = getopt(varargin, {'shallow_links'});
else
    args.shallow_links = 0; % Follow links
end

if ~exist(filename)
    error('File does not exist at the specified path')
end

% Open the associated file for writing
fp = fopen(filename, 'r');

header = fhd_read_header(fp);
data = [];
data = fhd_read_group(fp, data, args);
data.header = header;

fclose(fp);

end

%element_type_group = 0;
%element_type_attribute = 1;
%element_type_dataset = 2;

function [header] = fhd_read_header(fp)
% Read the header from the specified file pointer
    temp = fread(fp, 3, 'uint8=>char');
    if strcmp(char(reshape(temp, 1, 3)), 'fhd') ~= 1
        error('FHD header is invalid')
    end
    header = [];
    fread(fp, 1, 'uint8'); % Read null character
    header.major_version = fread(fp, 1, 'uint8');
    header.minor_version = fread(fp, 1, 'uint8');
    header.minor_minor_version = fread(fp, 1, 'uint8');
    header.pointer_size = fread(fp, 1, 'uint8');
    header.num_pointer_entries = fread(fp, 1, 'uint64');
    
    % Seek to end of header
    fseek(fp, 1024, 'bof');
end

function [data] = fhd_read_group(fp, data, args)
    current_location = ftell(fp);
    
    group_location = fread(fp, 1, 'uint64');
    if group_location ~= current_location
        error('Group location is invalid')
    end
    element_type = fread(fp, 1, 'uint8');
    if (element_type ~= 0)
        error('Element type is invalid for group')
    end
    name_length = fread(fp, 1, 'uint16');
    name = fread(fp, name_length, 'uint8=>char');

    name = char(reshape(name, 1, name_length));
    parent = fread(fp, 1, 'uint64');
    
    % Load children
    if name == '/'
        data = fhd_read_linked_list(fp, data, args);
    else
        data.(name) = [];
        data.(name) = fhd_read_linked_list(fp, [], args);
    end
end

function [data] = fhd_read_linked_list(fp, data, args)
    next_location = fread(fp, 1, 'uint64');
    previous_location = fread(fp, 1, 'uint64');
    max_entries = fread(fp, 1, 'uint64');
    num_contents = fread(fp, 1, 'uint64');
    if num_contents > max_entries
        error('Number of contents for linked list exceeds max entries')
    end
    child_locations = fread(fp, num_contents, 'uint64');
    
    % Load each of our children
    for i = 1:length(child_locations)
        fseek(fp, child_locations(i), 'bof');
        child_location = fread(fp, 1, 'uint64');
        if child_location ~= child_locations(i)
            error('Child location is invalid')
        end
        child_element_type = fread(fp, 1, 'uint8');
        fseek(fp, child_locations(i), 'bof');
        if child_element_type == 0
            data = fhd_read_group(fp, data, args);
        elseif child_element_type == 1
            data = fhd_read_attribute(fp, data, args);
        else
            data = fhd_read_dataset(fp, data, args);
        end
    end
    
    % Load the next set of entries in our list
    if (next_location > 0)
        fseek(fp, next_location, 'bof');
        data = fhd_read_linked_list(fp, data, args);
    end
end

function [type] = fhd_datatype(data_type)
    if (data_type == 0)
        type = 'double';
    elseif (data_type == 1)
        type = 'single';
    elseif (data_type == 2)
        type = 'uint8';
    elseif (data_type == 3)
        type = 'uint16';
    elseif (data_type == 4)
        type = 'uint32';
    elseif (data_type == 5)
        type = 'uint64';
    elseif (data_type == 6)
        type = 'int8';
    elseif (data_type == 7)
        type = 'int16';
    elseif (data_type == 8)
        type = 'int32';
    elseif (data_type == 9)
        type = 'int64';
    elseif (data_type == 10)
        type = 'uint8=>char';
    elseif (data_type == 11)
        type = 'uint64';
    else
        error('Unknown datatype')
    end
end

function [data] = fhd_read_attribute(fp, data, args)
    current_location = ftell(fp);
    location = fread(fp, 1, 'uint64');
    if location ~= current_location
        error('Attribute location is invalid')
    end
    element_type = fread(fp, 1, 'uint8');
    if (element_type ~= 1)
        error('Element type is invalid for attribute')
    end
    name_length = fread(fp, 1, 'uint16');
    name = fread(fp, name_length, 'uint8=>char');
    name = char(reshape(name, 1, name_length));
    parent = fread(fp, 1, 'uint64');
    num_dimensions = fread(fp, 1, 'uint8');
    dimensions = fread(fp, num_dimensions, 'uint64');
    data_type = fread(fp, 1, 'uint8');
    temp = fread(fp, prod(dimensions), fhd_datatype(data_type));
    if num_dimensions > 1
        temp = reshape(data, dimensions(:));
    end
    if data_type == 10
        temp = reshape(temp, 1, numel(temp));
    end
    
    data.(name) = temp;
end

function [data] = fhd_read_dataset(fp, data, args)
    current_location = ftell(fp);
    location = fread(fp, 1, 'uint64');
    if location ~= current_location
        error('Dataset location is invalid')
    end
    element_type = fread(fp, 1, 'uint8');
    if (element_type ~= 2)
        error('Element type is invalid for dataset')
    end
    name_length = fread(fp, 1, 'uint16');
    name = fread(fp, name_length, 'uint8=>char');
    name = char(reshape(name, 1, name_length));
    parent = fread(fp, 1, 'uint64');
    num_dimensions = fread(fp, 1, 'uint8');
    dimensions = fread(fp, num_dimensions, 'uint64');
    data_type = fread(fp, 1, 'uint8');
    long_dimension = fread(fp, 1, 'uint64');
    data.(name) = zeros(dimensions(:), long_dimension);
    data.(name) = fhd_read_dataset_linked_list(fp, dimensions, long_dimension, data_type, args);
end

function [data] = fhd_read_dataset_linked_list(fp, dimensions, long_dimension, data_type, args)
    if data_type ~= 11 || args.shallow_links % Is not a pointer
        data = zeros(dimensions(:), long_dimension, fhd_datatype(data_type));
    else
        data = cell(dimensions(:), long_dimension);
    end
    index = 1;
    while (1)
        next_location = fread(fp, 1, 'uint64');
        previous_location = fread(fp, 1, 'uint64');
        max_entries = fread(fp, 1, 'uint64');
        num_contents = fread(fp, 1, 'uint64');
        if num_contents > max_entries
            error('Number of contents for linked list exceeds max entries')
        end
        child_locations = fread(fp, num_contents, 'uint64');

        % Load each of our children
        for i = 1:length(child_locations)
            fseek(fp, child_locations(i), 'bof');
            long_axis = fread(fp, 1, 'uint64');
            temp = fread(fp, prod(dimensions) * long_axis, fhd_datatype(data_type));
            otherdims = repmat({':'}, 1, length(dimensions));
            if (data_type == 11 && ~args.shallow_links)
                temp_pointers = cell(1, numel(temp));
                for j = 1:length(temp_pointers)
                    fseek(fp, temp(j), 'bof');
                    current_group = fhd_read_group(fp, [], args);
                    fields = fieldnames(current_group);
                    temp_pointers{j} = current_group.(fields{1});
                    temp_pointers{j}.name = fields{1};
                end
                temp = temp_pointers;
            end
            temp = reshape(temp, dimensions(:), long_axis);
            data(otherdims{:}, index:index+long_axis-1) = temp;
            index = index + long_axis;
        end
        if next_location == 0
            break;
        else
            fseek(fp, next_location, 'bof');
        end
    end
end

