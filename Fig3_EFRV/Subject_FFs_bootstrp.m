%% Subject FFs: TEST
%%%%%%%
% R29 %
%%%%%%%
disp('R29_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
[file_list1] = find_files('R29','17_32');

R29_FFs = cell(length(file_list1),1);
R29_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [R29_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    R29_N(p)              = length(R29_FFs{p});
end

% Determine minimum number of units for subject
R29_N   = min(R29_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(R29_FFs)
    
    data   = R29_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),R29_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

R29_mean = mean(bs_mean);
R29_SE   = std(bs_mean);

%%%%%%%
% V1 %
%%%%%%%
disp('V1_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
[file_list2] = find_files('V1','1_16');

V1_FFs = cell(length(file_list2),1);
V1_N   = zeros(length(file_list2),1);

for p = 1:length(file_list2)
    file = file_list2{p};
    [V1_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V1_N(p)              = length(V1_FFs{p});
end

% Determine minimum number of units for subject
V1_N   = min(V1_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list2));

for p = 1:length(V1_FFs)
    
    data   = V1_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V1_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V1_mean = mean(bs_mean);
V1_SE   = std(bs_mean);

%%%%%%%
% V4 %
%%%%%%%
disp('V4_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
[file_list3] = find_files('V4','1_16');

V4_FFs = cell(length(file_list3),1);
V4_N   = zeros(length(file_list3),1);

for p = 1:length(file_list3)
    file = file_list3{p};
    [V4_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V4_N(p)              = length(V4_FFs{p});
end

% Determine minimum number of units for subject
V4_N   = min(V4_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list3));

for p = 1:length(V4_FFs)
    
    data   = V4_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V4_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V4_mean = mean(bs_mean);
V4_SE   = std(bs_mean);

%%%%%%%
% V7 %
%%%%%%%
disp('V7_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
[file_list4] = find_files('V7','1_16');

V7_FFs = cell(length(file_list4),1);
V7_N   = zeros(length(file_list4),1);

for p = 1:length(file_list4)
    file = file_list4{p};
    [V7_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V7_N(p)              = length(V7_FFs{p});
end

% Determine minimum number of units for subject
V7_N   = min(V7_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list4));

for p = 1:length(V7_FFs)
    
    data   = V7_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V7_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V7_mean = mean(bs_mean);
V7_SE   = std(bs_mean);

%%%%%%%
% V8 %
%%%%%%%
disp('V8_test')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list5] = find_files('V8','1_16');

V8_FFs = cell(length(file_list5),1);
V8_N   = zeros(length(file_list5),1);

for p = 1:length(file_list5)
    file = file_list5{p};
    [V8_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V8_N(p)              = length(V8_FFs{p});
end

% Determine minimum number of units for subject
V8_N   = min(V8_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list5));

for p = 1:length(V8_FFs)
    
    data   = V8_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V8_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V8_mean = mean(bs_mean);
V8_SE   = std(bs_mean);

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig3_EFRV')
save('Subject_FFS_test_bootstrp.mat','*_N','*_SE','*_mean','file_*')
clear

%% Subject FFs: CONTROL
%%%%%%%
% R29 %
%%%%%%%
disp('R29_control')
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
[file_list1] = find_files('R29','1_16');

R29_FFs = cell(length(file_list1),1);
R29_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [R29_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    R29_N(p)              = length(R29_FFs{p});
end

% Determine minimum number of units for subject
R29_N   = min(R29_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(R29_FFs)
    
    data   = R29_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),R29_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

R29_mean = mean(bs_mean);
R29_SE   = std(bs_mean);

%%%%%%%
% V1 %
%%%%%%%
disp('V1_control')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
[file_list2] = find_files('V1','17_32');

V1_FFs = cell(length(file_list2),1);
V1_N   = zeros(length(file_list2),1);

for p = 1:length(file_list2)
    file = file_list2{p};
    [V1_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V1_N(p)              = length(V1_FFs{p});
end

% Determine minimum number of units for subject
V1_N   = min(V1_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list2));

for p = 1:length(V1_FFs)
    
    data   = V1_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V1_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V1_mean = mean(bs_mean);
V1_SE   = std(bs_mean);

%%%%%%%
% V4 %
%%%%%%%
disp('V4_control')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
[file_list3] = find_files('V4','17_32');

V4_FFs = cell(length(file_list3),1);
V4_N   = zeros(length(file_list3),1);

for p = 1:length(file_list3)
    file = file_list3{p};
    [V4_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V4_N(p)              = length(V4_FFs{p});
end

% Determine minimum number of units for subject
V4_N   = min(V4_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list3));

for p = 1:length(V4_FFs)
    
    data   = V4_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V4_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V4_mean = mean(bs_mean);
V4_SE   = std(bs_mean);


%%%%%%%
% V7 %
%%%%%%%
disp('V7_control')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
[file_list4] = find_files('V7','17_32');

V7_FFs = cell(length(file_list4),1);
V7_N   = zeros(length(file_list4),1);

for p = 1:length(file_list4)
    file = file_list4{p};
    [V7_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V7_N(p)              = length(V7_FFs{p});
end

% Determine minimum number of units for subject
V7_N   = min(V7_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list4));

for p = 1:length(V7_FFs)
    
    data   = V7_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V7_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V7_mean = mean(bs_mean);
V7_SE   = std(bs_mean);

%%%%%%%
% V8 %
%%%%%%%
disp('V8_control')
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list5] = find_files('V8','17_32');

V8_FFs = cell(length(file_list5),1);
V8_N   = zeros(length(file_list5),1);

for p = 1:length(file_list5)
    file = file_list5{p};
    [V8_FFs{p},Behavior] = Tracking_FFs(file,1000,'test');
    V8_N(p)              = length(V8_FFs{p});
end

% Determine minimum number of units for subject
V8_N   = min(V8_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list5));

for p = 1:length(V8_FFs)
    
    data   = V8_FFs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V8_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V8_mean = mean(bs_mean);
V8_SE   = std(bs_mean);

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig3_EFRV')
save('Subject_FFs_control_bootstrp.mat','*_N','*_SE','*_mean','file_*')
% Shut the system down when done computing
system('%windir%\system32\rundll32.exe PowrProf.dll, SetSuspendState')
