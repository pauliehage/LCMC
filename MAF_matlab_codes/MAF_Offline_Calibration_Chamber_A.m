clc;
clear;

file_path = [pwd filesep];
[file_name_cal, file_path_cal] = uigetfile([file_path '*.JSON'],...
    'Select calibration.JSON file');

fid = fopen([file_path_cal, file_name_cal]); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 

[file_name_C, file_path_C] = uigetfile([file_path '*.mat'],...
    'Select ANALYZED_RECAL.mat file');

load([file_path_C, file_name_C],'CAL_MATRIX');
C = CAL_MATRIX;


value = jsondecode(str);

P = [value.right_calibration.p0,...
    value.right_calibration.p1,...
    value.right_calibration.p2;...
    value.right_calibration.p3,...
    value.right_calibration.p4,...
    value.right_calibration.p5];

Q = P;

Q(1:2,1:2) = C(1:2,1:2)'*P(1:2,1:2);
Q(:,3) = C(:,1:2)' * [P(:,3);1];

value.right_calibration.p0 = Q(1,1);
value.right_calibration.p1 = Q(1,2);
value.right_calibration.p2 = Q(1,3);
value.right_calibration.p3 = Q(2,1);
value.right_calibration.p4 = Q(2,2);
value.right_calibration.p5 = Q(2,3);

value_recal = jsonencode(value,'PrettyPrint',true);
fid = fopen([file_path_cal, 'recal.JSON'],'W'); 
fwrite(fid,value_recal);
fclose(fid); 