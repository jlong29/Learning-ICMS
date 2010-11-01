%% Subject SDs: TEST
%%%%%%%
% R29 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
[file_list1] = find_files('R29','17_32');

R29_SDs = cell(length(file_list1),1);
R29_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [R29_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    R29_N(p)              = length(R29_SDs{p});
end

% Determine minimum number of units for subject
R29_N   = min(R29_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(R29_SDs)
    
    data   = R29_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),R29_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

R29_mean = mean(bs_mean);
R29_SE   = std(bs_mean);

figure
bar(R29_mean,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       errorbar(1:length(file_list1),R29_mean,-R29_SE,R29_SE,'color','k','linestyle','none','linewidth',3)
       title('R29')

%%%%%%%
% V1 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
[file_list2] = find_files('V1','1_16');

V1_SDs = cell(length(file_list1),1);
V1_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [V1_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V1_N(p)              = length(V1_SDs{p});
end

% Determine minimum number of units for subject
V1_N   = min(V1_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(V1_SDs)
    
    data   = V1_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V1_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V1_mean = mean(bs_mean);
V1_SE   = std(bs_mean);

figure
bar(V1_mean,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       errorbar(1:length(file_list2),V1_mean,-V1_SE,V1_SE,'color','k','linestyle','none','linewidth',3)
       title('V1')

%%%%%%%
% V4 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
[file_list3] = find_files('V4','1_16');

V4_SDs = cell(length(file_list1),1);
V4_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [V4_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V4_N(p)              = length(V4_SDs{p});
end

% Determine minimum number of units for subject
V4_N   = min(V4_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(V4_SDs)
    
    data   = V4_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V4_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V4_mean = mean(bs_mean);
V4_SE   = std(bs_mean);

figure
bar(V4_mean,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       errorbar(1:length(file_list3),V4_mean,-V4_SE,V4_SE,'color','k','linestyle','none','linewidth',3)
       title('V4')

%%%%%%%
% V7 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
[file_list4] = find_files('V7','1_16');

V7_SDs = cell(length(file_list1),1);
V7_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [V7_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V7_N(p)              = length(V7_SDs{p});
end

% Determine minimum number of units for subject
V7_N   = min(V7_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(V7_SDs)
    
    data   = V7_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V7_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V7_mean = mean(bs_mean);
V7_SE   = std(bs_mean);

figure
bar(V7_mean,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       errorbar(1:length(file_list4),V7_mean,-V7_SE,V7_SE,'color','k','linestyle','none','linewidth',3)
       title('V7')

%%%%%%%
% V8 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list5] = find_files('V8','1_16');

V8_SDs = cell(length(file_list1),1);
V8_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [V8_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V8_N(p)              = length(V8_SDs{p});
end

% Determine minimum number of units for subject
V8_N   = min(V8_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(V8_SDs)
    
    data   = V8_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V8_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V8_mean = mean(bs_mean);
V8_SE   = std(bs_mean);

figure
bar(V8_mean,'facecolor',[0.6902 0.7137 0.6118],'edgecolor','k'),hold on,
       errorbar(1:length(file_list5),V8_mean,-V8_SE,V8_SE,'color','k','linestyle','none','linewidth',3)
       title('V8')

%% Subject SDs: CONTROL
%%%%%%%
% R29 %
%%%%%%%
cd('C:\Data\Plexon_Data\Test_Data\Learning\R29')
[file_list1] = find_files('R29','1_16');

R29_SDs = cell(length(file_list1),1);
R29_N   = zeros(length(file_list1),1);

for p = 1:length(file_list1)
    file = file_list1{p};
    [R29_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    R29_N(p)              = length(R29_SDs{p});
end

% Determine minimum number of units for subject
R29_N   = min(R29_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list1));

for p = 1:length(R29_SDs)
    
    data   = R29_SDs{p};
    
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
cd('C:\Data\Plexon_Data\Test_Data\Learning\V1')
[file_list2] = find_files('V1','17_32');

V1_SDs = cell(length(file_list2),1);
V1_N   = zeros(length(file_list2),1);

for p = 1:length(file_list2)
    file = file_list2{p};
    [V1_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V1_N(p)              = length(V1_SDs{p});
end

% Determine minimum number of units for subject
V1_N   = min(V1_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list2));

for p = 1:length(V1_SDs)
    
    data   = V1_SDs{p};
    
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
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4')
[file_list3] = find_files('V4','17_32');

V4_SDs = cell(length(file_list3),1);
V4_N   = zeros(length(file_list3),1);

for p = 1:length(file_list3)
    file = file_list3{p};
    [V4_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V4_N(p)              = length(V4_SDs{p});
end

% Determine minimum number of units for subject
V4_N   = min(V4_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list3));

for p = 1:length(V4_SDs)
    
    data   = V4_SDs{p};
    
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
cd('C:\Data\Plexon_Data\Test_Data\Learning\V7')
[file_list4] = find_files('V7','17_32');

V7_SDs = cell(length(file_list4),1);
V7_N   = zeros(length(file_list4),1);

for p = 1:length(file_list4)
    file = file_list4{p};
    [V7_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V7_N(p)              = length(V7_SDs{p});
end

% Determine minimum number of units for subject
V7_N   = min(V7_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list4));

for p = 1:length(V7_SDs)
    
    data   = V7_SDs{p};
    
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
cd('C:\Data\Plexon_Data\Test_Data\Learning\V8')
[file_list5] = find_files('V8','17_32');

V8_SDs = cell(length(file_list5),1);
V8_N   = zeros(length(file_list5),1);

for p = 1:length(file_list5)
    file = file_list5{p};
    [V8_SDs{p},Behavior] = Tracking_SDs(file,1000,'test');
    V8_N(p)              = length(V8_SDs{p});
end

% Determine minimum number of units for subject
V8_N   = min(V8_N);
bs_N    = 500;
bs_mean = zeros(bs_N,length(file_list5));

for p = 1:length(V8_SDs)
    
    data   = V8_SDs{p};
    
    for q = 1:bs_N
        rs         = randsample(length(data),V8_N,'true');
    
        bs_mean(q,p) = mean(data(rs));
    
    end
    
end

V8_mean = mean(bs_mean);
V8_SE   = std(bs_mean);

