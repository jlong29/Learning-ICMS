% Reviewer#1 brings up the valid concern that the observed changes in ensemble
% firing rate variability may be attributable to changes in the ensemble firing
% rate between sessions. This is reasonable given the roughly poisson nature of
% the data, which means the variance is equal to the meaning firing rate.
% Therefore, an increase the unit firing rate would result in an increase
% variance.
close all
set(0,'units','inches')
scrnsz = get(0,'Screensize');
close all
figure('units','inches','Position',[1 scrnsz(4)*(2/3) 8.5 4]) 

% basic relationship of among the statistics
lambda = 3;
win   = 100;
x     = poissrnd(lambda,100*win,1);
F     = zeros(100,1);
u     = zeros(100,1);
sig2  = zeros(100,1);

for n = 1:length(x)/win
    u(n)    = mean(x(1 + win*(n-1): win*n));
    sig2(n) = var(x(1 + win*(n-1): win*n));
    F(n)    = sig2(n)/u(n);
end

subplot(121), plot(F,'b','linewidth',2),hold on, plot(u,'r','linewidth',2),...
    plot(sig2,'g','linewidth',2)
line([0 length(F)],[1 1],'linestyle','--','color','k')
line([0 length(F)],[lambda lambda],'linestyle','--','color','k')

xlabel(['Repititions of ' num2str(win) ' samples'])
ylabel('Real Number Values')
title('Behavior of \mu, \sigma^2, and Fano Factor')

ylim([0 max([u;sig2;F])+.5*lambda])
text(5,max([u;sig2;F])+.4*lambda,...
    [{['lambda = ' num2str(lambda)]};...
    {['win = ' num2str(win)]};...
    {['Corr(\mu,\sigma^2) = ' num2str(corr(u,sig2))]};...
    {['Corr(\mu,F) = ' num2str(corr(u,F))]};...
    {['Corr(F,\sigma^2) = ' num2str(corr(F,sig2))]}],'VerticalAlignment','top')

legend('Fano Factor','\mu','\sigma^2')

% The basic phenomenon and how
N      = 100;
lambda = 10;
x      = poissrnd(lambda,1,N);

y      = poissrnd(lambda + 10,1,N);

subplot(122),plot(x,'k','linewidth',2),hold on, plot(y,'linewidth',2,'color',[0.624   0.624   0.624])
xlabel('Samples')
ylabel('Counts')
title([{'Relationship between mean, variance,'};{'and Fano Factor'}])

ylim([0 max([x y])+.5*lambda])

F1    = var(x)/mean(x);
F2    = var(y)/mean(y);

text(5,max([x y])+.4*lambda,...
    [{['\mu_1 = ' num2str(mean(x)) '  \sigma_1^2 = ' num2str(var(x))]};...
     {['\mu_2 = ' num2str(mean(y)) '  \sigma_2^2 = ' num2str(var(y))]};...
     {['F_1 = ' num2str(F1)]};...
     {['F_2 = ' num2str(F2)]}],'VerticalAlignment','top')

legend('data_1','data_2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Global figure Settings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(findobj(gcf,'Type','axes'),'box','off')
set(gcf,'color','w')

set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 8.5 4])
