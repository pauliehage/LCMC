function PGH_csv_to_mat(folder_path)
% get the folder path and convert all the csv files in that folder to mat files

% if there is no inputs, then set folder_path to pwd
if nargin < 1
    folder_path = pwd;
end
% add filesep to the end of folder_path
if ~strcmp(folder_path(end), filesep)
    folder_path = [folder_path filesep];
end
% get the list of fhd files to analyze
FILES_CSV = dir([folder_path '*.csv']);

for counter_files = 1 : length(FILES_CSV)
    file_csv = FILES_CSV(counter_files);
    [~,file_name,~] = fileparts(file_csv.name);
    file_path_name = [folder_path file_csv.name];

    fprintf([num2str(counter_files) '. Loading: ' file_name ' ... '])

    clearvars('data');

    data=csv_load(file_path_name);

    save([folder_path filesep file_name '.mat'], 'data' );

    fprintf(' --> Completed.\n')
end
fprintf('ALL Completed. \n')
end

function [data] = csv_load(filename)

if ~exist(filename)
    error('File does not exist at the specified path')
end

% load raw table from csv file
csv_table = xlsread(filename);

% load columns into variable
tip_tongue_x = csv_table(:,2);
tip_tongue_y = csv_table(:,3);
r_tongue_x = csv_table(:,5);
r_tongue_y = csv_table(:,6);
l_tongue_x = csv_table(:,8);
l_tongue_y = csv_table(:,9);
mid_tongue_x = csv_table(:,11);
mid_tongue_y = csv_table(:,12);
r_nose_x = csv_table(:,14);
r_nose_y = csv_table(:,15);
l_nose_x = csv_table(:,17);
l_nose_y = csv_table(:,18);
r_food_x = csv_table(:,20);
r_food_y = csv_table(:,21);
l_food_x = csv_table(:,23);
l_food_y = csv_table(:,24);
r_tube_r_x = csv_table(:,26);
r_tube_r_y = csv_table(:,27);
r_tube_l_x = csv_table(:,29);
r_tube_l_y = csv_table(:,30);
l_tube_r_x = csv_table(:,32);
l_tube_r_y = csv_table(:,33);
l_tube_l_x = csv_table(:,35);
l_tube_l_y = csv_table(:,36);

% build table
data = table(tip_tongue_x,tip_tongue_y,r_tongue_x, r_tongue_y, ...
    l_tongue_x, l_tongue_y, mid_tongue_x, mid_tongue_y, r_nose_x, r_nose_y, ...
    l_nose_x, l_nose_y, r_food_x, r_food_y, l_food_x, l_food_y, r_tube_r_x, r_tube_r_y, ...
    r_tube_l_x, r_tube_l_y, l_tube_r_x, l_tube_r_y, l_tube_l_x, l_tube_l_y);

end