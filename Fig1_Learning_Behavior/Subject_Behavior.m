cd('C:\Data\Plexon_Data\Test_Data\Learning\Tracking_KLs_test\R29')
[file_list1] = find_files('R29','17_32_justICMS');

cd('C:\Data\Plexon_Data\Test_Data\Learning\Tracking_KLs_test\V1')
[file_list2] = find_files('V1','_1_16_full');

cd('C:\Data\Plexon_Data\Test_Data\Learning\Tracking_KLs_test\V7')
[file_list3] = find_files('V7','1_16_full');

cd('C:\Data\Plexon_Data\Test_Data\Learning\Tracking_KLs_test\V8')
[file_list4] = find_files('V8','1_16_full');

file_lists    = cell(4,1);
file_lists{1} = file_list1;
file_lists{2} = file_list2;
file_lists{3} = file_list3;
file_lists{4} = file_list4;

% Initialize all Subjects' Behavior variable
subject_Behavior = cell(5,1);

% for subjects: R29,V1,V7,V8
for p = 1:length(file_lists)
    file_list = file_lists{p};
    N         = length(file_list);

    for n = 1:N
        load(file_list{n})
        eval(['Behavior' num2str(n) ' = Behavior;'])
        eval(['behavior_ts' num2str(n) ' = behavior_ts;'])
    end

    % Finish: Determine range for KL-div over days %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    beh_sum = zeros(3,N);
    
    for n = 1:N
        % Grab current variables
        eval(['Behavior    = Behavior' num2str(n) ';'])
        
        %Calculate percentages
        bsums        = sum(diff(Behavior,[],2),2);
        total        = sum(bsums);
        
        beh_sum(:,n) = 100.*bsums./total;
    end
    subject_Behavior{p} = beh_sum;
end

% For subject V4
cd('C:\Data\Plexon_Data\Test_Data\Learning\V4\')
[file_list] = find_files('V4','_1_16');

beh_sum = zeros(3,length(file_list));
for n = 1:length(file_list)
    load(file_list{n})
    event = structpack(who('Eve*'));
    
    [event,Behavior,behavior_ts,Response_t] = Behavior_preproc(event);
    
    %Calculate percentages
    bsums        = sum(diff(Behavior,[],2),2);
    total        = sum(bsums);

    beh_sum(:,n) = 100.*bsums./total;
end

subject_Behavior{5} = beh_sum;
% Reorder
subject_Behavior = subject_Behavior([1 2 5 3 4]);
keep subject_Behavior