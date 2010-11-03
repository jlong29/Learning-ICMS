%%% Figure 7: SNR models to distance
close all
clear

%%%%%%%%%%%%%%%%%%%%%
%%% Figure layout %%%
%%%%%%%%%%%%%%%%%%%%%
% Set dots per inch
dpi = 96;
figure('units','pixels','Position',[1 1 8.5 7.6].*dpi)

% From left to right
H = [1 2.9 4.8];

% From bottom to top
V = [1.7 4.0];

% Layout of axes
pos      = [H(1)  V(2)  1.75 1.75
            H(2)  V(2)  1.75 1.75
            H(3)  V(2)  1.75 1.75
            
            H(1)  V(1)  1.75 1.75
            H(2)  V(1)  1.75 1.75
            H(3)  V(1)  1.75 1.75].*dpi;



% Data Directory        
cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig7_SNR_distance')
% Subject IDs
subjects = {'R29*','V1*','V4*','V7*','V8*'};
load('Fig7_SNR_distance_test_sub2')

% Grab the data
vars1 = who(subjects{1});
vars2 = who(subjects{2});
vars3 = who(subjects{3});
vars4 = who(subjects{4});
vars5 = who(subjects{5});
    
% Grab max and min for signal and internal noise
maxs  = zeros(length(vars1),2);
mins  = zeros(length(vars1),2);

for n = 1:length(vars1)
    temp1    = eval(vars1{n});
    temp2    = eval(vars2{n});
    temp3    = eval(vars3{n});
    temp4    = eval(vars4{n});
    temp5    = eval(vars5{n});
    
    temp1(isinf(temp1)) = NaN;
    temp2(isinf(temp2)) = NaN;
    temp3(isinf(temp3)) = NaN;
    temp4(isinf(temp4)) = NaN;
    temp5(isinf(temp5)) = NaN;
    
    maxs(n,:) = nanmax([temp1(:,1:2);temp2(:,1:2);temp3(:,1:2);temp4(:,1:2);temp5(:,1:2)]);
    mins(n,:) = nanmin([temp1(:,1:2);temp2(:,1:2);temp3(:,1:2);temp4(:,1:2);temp5(:,1:2)]);
end

% Calculate the max and min for each metric
ymin = floor(min(mins));
ymax = ceil(max(maxs));

% temp to correct for outlier
ymax(1) = 50;

for n = 1:length(vars1)
    
    % Data for each subject
    temp1 = eval(vars1{n});
    temp2 = eval(vars2{n});
    temp3 = eval(vars3{n});
    temp4 = eval(vars4{n});
    temp5 = eval(vars5{n});
    
    ind1  = ~logical(temp1(:,4));
    ind2  = logical(temp2(:,4));
    ind3  = logical(temp3(:,4));
    ind4  = logical(temp4(:,4));
    ind5  = logical(temp5(:,4));
    
    % Test Hemisphere data for each subject
    sub1  = temp1(ind1,:);
    [a,b] = sort(sub1(:,3),1);
    sub1  = sub1(b,:);
    
    sub2  = temp2(ind2,:);
    [a,b] = sort(sub2(:,3),1);
    sub2  = sub2(b,:);
    
    sub3  = temp3(ind3,:);
    [a,b] = sort(sub3(:,3),1);
    sub3  = sub3(b,:);
    
    sub4  = temp4(ind4,:);
    [a,b] = sort(sub4(:,3),1);
    sub4  = sub4(b,:);
    
    sub5  = temp5(ind5,:);
    [a,b] = sort(sub5(:,3),1);
    sub5  = sub5(b,:);
    
    h1 = axes('units','pixels','position',pos(n,:));
        plot(sub1(:,3),sub1(:,1),'+k'),hold on,
        plot(sub2(:,3),sub2(:,1),'ok'),
        plot(sub3(:,3),sub3(:,1),'sk'),
        plot(sub4(:,3),sub4(:,1),'vk'),
        plot(sub5(:,3),sub5(:,1),'*k'),
        
        % fit a line to the group data
        h = get(gca,'children');
        x = cell2mat(get(h,'xdata')');
        y = cell2mat(get(h,'ydata')');
        p = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        x = 0:.2:2.25;
        y = polyval(p,x);
        plot(x,y,'--k','linewidth',2)
        
        ylim([ymin(1) ymax(1)])
        xlim([0 2.25])
        axis square
        
        if n == 1
            title([{'Stage1:'};{'ICMS Signal'}],'fontname','arial','fontsize',10)
            ylabel('z-score','fontname','arial','fontsize',10)
            legend('S1','S2','S3','S4','S5','location','northwest')
        elseif n == 2
            title([{'Stage2:'};{'ICMS Signal'}],'fontname','arial','fontsize',10)
            set(gca,'yticklabel','')
        else
            title([{'Stage3:'};{'ICMS Signal'}],'fontname','arial','fontsize',10)
            set(gca,'yticklabel','')
        end
        
    h2 = axes('units','pixels','position',pos(n+3,:));
        plot(sub1(:,3),sub1(:,2),'+k'),hold on,
        plot(sub2(:,3),sub2(:,2),'ok'),
        plot(sub3(:,3),sub3(:,2),'sk'),
        plot(sub4(:,3),sub4(:,2),'vk'),
        plot(sub5(:,3),sub5(:,2),'*k'),
        
        % fit a line to the group data
        h = get(gca,'children');
        x = cell2mat(get(h,'xdata')');
        y = cell2mat(get(h,'ydata')');
        p = polyfit(x(~isnan(y)),y(~isnan(y)),1);
        x = 0:.2:2.25;
        y = polyval(p,x);
        plot(x,y,'--k','linewidth',2)
        
        ylim([ymin(2) ymax(2)])
        xlim([0 2.25])
        title(['Day ' num2str(n) ': internal noise'])
        xlabel([{'Distance from'}; {'stimulating electrode (mm)'}],'fontname','arial','fontsize',10)
        axis square

        if n == 1
            title('ICMS Internal Noise','fontname','arial','fontsize',10)
            ylabel('t-value','fontname','arial','fontsize',10)
        elseif n == 2
            title('ICMS Internal Noise','fontname','arial','fontsize',10)
            set(gca,'yticklabel','')
        else
            title('ICMS Internal Noise','fontname','arial','fontsize',10)
            set(gca,'yticklabel','')
        end
        
end

set(gcf,'color','w')
set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 7.6])
