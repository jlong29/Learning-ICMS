close all
clear

% Load Ensemble Firing Rate Variability Data
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
load('Subject_SDs_test.mat')

% Grab all per subject Accelerometer files
cd('C:\Data\Plexon_Data\Test_Data\R29\R29_Accelerometer')
[file_Acc1] = find_files('R29','_Acc');

cd('C:\Data\Plexon_Data\Test_Data\V1\V1_Accelerometer')
[file_Acc2] = find_files('V1','_Acc');

cd('C:\Data\Plexon_Data\Test_Data\V4\V4_Accelerometer')
[file_Acc3] = find_files('V4','_Acc');

cd('C:\Data\Plexon_Data\Test_Data\V7\V7_Accelerometer')
[file_Acc4] = find_files('V7','_Acc');

cd('C:\Data\Plexon_Data\Test_Data\V8\V8_Accelerometer')
[file_Acc5] = find_files('V8','_Acc');

subject_Acc = cell(5,2);

for n = 1:5
    % Grab files for subject
    files     = eval(['file_Acc' num2str(n)]);
    % Find subject ID
    [a,temp]  = fileparts(files{1});
    chk       = strfind(temp,'_');
    sub_ID    = temp(1:chk(1)-1);
    
    % Generate Acceleration Measure and standard deviation
    normA_std = eval(['zeros(length(file_Acc' num2str(n) '),1)']);

    for p = 1:length(files)
        load(files{p})
        Acc = who('AD*');

        temp = zeros(length(eval(Acc{1})),length(Acc));
        for q = 1:length(Acc)
            temp(:,q) = abs(zscore(eval(Acc{q})));
        end

        eval(['normA' num2str(p) '= mean(temp,2);'])

        normA_std(p) = eval(['std(normA' num2str(p) ');']);
        clear AD* a*
    end

    % Calculate Threshold Events
    thrs   = zeros(6,1);
    thrsev = zeros(length(thrs),length(files));

    for p = 1:length(files);

        thrs = normA_std(p)*5:normA_std(p):normA_std(p)*10;

        for q = 1:length(thrs)
            temp        = eval(['normA' num2str(p)]);
            t           = find(temp > thrs(q));
            % Filter out consecutive time points
            t           = diff(t) > 1;
            thrsev(q,p) = sum(t);
        end
    end
    
    % Update per subject Acceleration event array
    subject_Acc{n,1} = thrsev;
    
    % Collate corresponding Ensemble Firing Rate Variability
    ind      = zeros(length(files),1);
        
    for p = 1:length(files)
        [a,file] = fileparts(files{p});
        file     = file(1:end-4);
        
        fileEfrV = eval(['file_list' num2str(n)]);
        for q = 1:length(fileEfrV)
            chk  = strfind(fileEfrV{q},file);
            if ~isempty(chk)
                ind(p) = q;
            end
        end
    end
    
    % Update per subject EfrV data
    Efrv             = eval([sub_ID '_mean;']);
    subject_Acc{n,2} = Efrv(ind(logical(ind)));
    
    % Adjust for differences in number of files
    if size(subject_Acc{n,1},2) ~= size(subject_Acc{n,2},2)
        subject_Acc{n,1} = subject_Acc{n,1}(:,ind(logical(ind)));
    end
end

cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
save('Null_Acc2EfrV_2.mat','subject_Acc')

%%
clear
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures')
load('Null_Acc2EfrV_2.mat')

for q = 1:6
    scrnsz = get(0,'screensize');
    figure('Position',scrnsz)
    
    % Each Subject with different Marker and Mark Active Exploration Sessions
    plot(subject_Acc{1,1}(q,:),subject_Acc{1,2},'ok','markerfacecolor','k','markersize',5),hold on
    plot(subject_Acc{2,1}(q,:),subject_Acc{2,2},'sk','markerfacecolor','k','markersize',5)
    plot(subject_Acc{3,1}(q,:),subject_Acc{3,2},'dk','markerfacecolor','k','markersize',5)
    plot(subject_Acc{4,1}(q,:),subject_Acc{4,2},'^k','markerfacecolor','k','markersize',5)
    plot(subject_Acc{5,1}(q,:),subject_Acc{5,2},'<k','markerfacecolor','k','markersize',5)
    
    plot(subject_Acc{1,1}(q,6),subject_Acc{1,2}(6),'or','markerfacecolor','r','markersize',5)
    plot(subject_Acc{2,1}(q,3),subject_Acc{2,2}(3),'sr','markerfacecolor','r','markersize',5)
    plot(subject_Acc{3,1}(q,4),subject_Acc{3,2}(4),'dr','markerfacecolor','r','markersize',5)
    plot(subject_Acc{4,1}(q,5),subject_Acc{4,2}(5),'^r','markerfacecolor','r','markersize',5)
    plot(subject_Acc{5,1}(q,4),subject_Acc{5,2}(4),'<r','markerfacecolor','r','markersize',5)
    legend('Subject1','Subject2','Subject3','Subject4','Subject5')
    axis square
    
    title(['Threshold = ' num2str(q)])
end