% Herein, I'm examining the variability of neural firing rates when the subjects
% were actively engaged in learning. These are the 'test' periods i.e. when only
% ICMS was presented. Not ICMS + Tone and not Tone alone.
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
[file_list1] = find_files('R29','17_32');

R29_SDs = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [ensemble_SDs,Behavior] = Tracking_SDs(file,1000,'test');
    R29_SDs(p) = mean(ensemble_SDs);
end

keep R29_SDs file_list1

cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
[file_list2] = find_files('V1','1_16');

V1_SDs = zeros(length(file_list2),1);

for p = 1:length(file_list2)
    file = file_list2{p};
    [ensemble_SDs,Behavior] = Tracking_SDs(file,1000,'test');
    V1_SDs(p) = mean(ensemble_SDs);
end

keep R29_SDs file_list1 V1_SDs file_list2

cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
[file_list3] = find_files('V4','1_16');

V4_SDs = zeros(length(file_list3),1);

for p = 1:length(file_list3)
    file = file_list3{p};
    [ensemble_SDs,Behavior] = Tracking_SDs(file,1000,'test');
    V4_SDs(p) = mean(ensemble_SDs);
end

keep R29_SDs file_list1 V1_SDs file_list2 V4_SDs file_list3

cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
[file_list4] = find_files('V7','1_16');

V7_SDs = zeros(length(file_list4),1);

for p = 1:length(file_list4)
    file = file_list4{p};
    [ensemble_SDs,Behavior] = Tracking_SDs(file,1000,'test');
    V7_SDs(p) = mean(ensemble_SDs);
end

keep R29_SDs file_list1 V1_SDs file_list2 V4_SDs file_list3 V7_SDs file_list4

cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list5] = find_files('V8','1_16');

V8_SDs = zeros(length(file_list5),1);

for p = 1:length(file_list5)
    file = file_list5{p};
    [ensemble_SDs,Behavior] = Tracking_SDs(file,1000,'test');
    V8_SDs(p) = mean(ensemble_SDs);
end

keep R29_SDs file_list1 V1_SDs file_list2 V4_SDs file_list3 V7_SDs file_list4 V8_SDs file_list5

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
save('Subject_SDs')

quit

