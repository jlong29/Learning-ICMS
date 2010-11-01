% Addressing the concern about ensemble firing rate varability being
% proportional to the length of the recording (IDEA DRAFT)
set(0,'units','inches')
scrnsz = get(0,'Screensize');
close all
figure('units','inches','Position',[1 scrnsz(4)*(2/3) 8.5 4]) 

% If the reviewer's concern is correct things should look like this...
Tsession = 15.*randn(40,1) + 45;
EFRV     = .1.*Tsession + 1.5.*randn(length(Tsession),1);

subplot(1,3,1),plot(Tsession, EFRV,'ok','markerfacecolor','k')
lin = polyfit(Tsession,EFRV,1);
hold on, plot(sort(Tsession),polyval(lin,sort(Tsession)),'--k')
text(min([Tsession;5]), max(EFRV)-1,['R = ' num2str(corr(Tsession,EFRV))]) 
title([{'Reviewer # 1, Comment 8:'};{'His prediction'}])
ylabel('Ensemble Firing Rate Variability')
xlabel('Session Length (mins)')
xlim([0 max(Tsession)+5])
axis square

% The cleanest difference would be...
EFRV    = 10.*rand(length(Tsession),1) + 2;
subplot(1,3,2),plot(Tsession, EFRV,'ok','markerfacecolor','k'),hold on
lin = polyfit(Tsession,EFRV,1);
plot(sort(Tsession),polyval(lin,sort(Tsession)),'--k')

text(min([Tsession;5]), max(EFRV)-1,['R = ' num2str(corr(Tsession,EFRV))]) 
title('Clearest Refutation of Concern')
xlabel('Session Length (mins)')
xlim([0 max(Tsession)+ 5])
axis square

% It is likely there is a slight linear trend in EFRV vs. Session length, but I
% anticipate that during exploration, this trend will not explain the increase
% in EFRV
temp     = mean(Tsession);
Tsession = Tsession(1:end-4);
EFRV     = .1.*Tsession + 1.5.*randn(length(Tsession),1);
tsession = temp+10+5.*randn(5,1);
eEFRV    = (max(EFRV)+4)+ randn(5,1);

subplot(1,3,3),plot(tsession,eEFRV,'oc','markerfacecolor','c'),hold on
plot(Tsession,EFRV,'ok','markerfacecolor','k')
Tsession = [Tsession; tsession];
EFRV     = [EFRV; eEFRV];
lin      = polyfit(Tsession,EFRV,1);
plot(sort(Tsession),polyval(lin,sort(Tsession)),'--k')

text(min([Tsession;5]), max(EFRV)-1,['R = ' num2str(corr(Tsession,EFRV))]) 
title([{'What I anticipate the data'};{'will look like. Cyan marks'};{'Exploratory Session'}])
ylabel('Ensemble Firing Rate Variability')
xlabel('Session Length (mins)')
xlim([0 max(Tsession)+ 5])
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Global figure Settings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')

set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 4])
