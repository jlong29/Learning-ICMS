cd('C:\Program Files\MATLAB\R2006b\work\Learning_ICMS_figures\Fig5_SNR_Models')
load('Learning_ICMS_SNR_models_new.mat')

a = histc(S(:,1),linspace(min(S(:)),max(S(:)),100));
b = histc(S(:,2),linspace(min(S(:)),max(S(:)),100));
c = histc(S(:,3),linspace(min(S(:)),max(S(:)),100));
figure,bar(a,'facecolor','b')
hold on, bar(b,'facecolor','r')
hold on, bar(c,'facecolor','g')
title('SNR Models: Signal')
legend('Stage1','Stage2','Stage3')

a = histc(NI(:,1),linspace(min(NI(:)),max(NI(:)),100));
b = histc(NI(:,2),linspace(min(NI(:)),max(NI(:)),100));
c = histc(NI(:,3),linspace(min(NI(:)),max(NI(:)),100));
figure,bar(a,'facecolor','b')
hold on, bar(b,'facecolor','r')
hold on, bar(c,'facecolor','g')
title('SNR Models: Internal Noise')
legend('Stage1','Stage2','Stage3')

% Test hemisphere
edge = linspace(min(NE(:)),max(NE(:)),100);
a = histc(NE(:,1),linspace(min(NE(:)),max(NE(:)),100));
b = histc(NE(:,2),linspace(min(NE(:)),max(NE(:)),100));
c = histc(NE(:,3),linspace(min(NE(:)),max(NE(:)),100));
figure,bar(edge,a,'facecolor','b')
hold on, bar(edge,b,'facecolor','r')
hold on, bar(edge,c,'facecolor','g')
title('SNR Models: External Noise Total')
legend('Stage1','Stage2','Stage3')

% edge = linspace(min(NE_e(:)),max(NE_e(:)),100);
% a = histc(NE_e(:,1),linspace(min(NE_e(:)),max(NE_e(:)),100));
% b = histc(NE_e(:,2),linspace(min(NE_e(:)),max(NE_e(:)),100));
% c = histc(NE_e(:,3),linspace(min(NE_e(:)),max(NE_e(:)),100));
% figure,bar(edge,a,'facecolor','b')
% hold on, bar(edge,b,'facecolor','r')
% hold on, bar(edge,c,'facecolor','g')
% title('SNR Models: External Noise Excitatory')
% legend('Stage1','Stage2','Stage3')
% 
% edge = linspace(min(NE_i(:)),max(NE_i(:)),100);
% a = histc(NE_i(:,1),linspace(min(NE_i(:)),max(NE_i(:)),100));
% b = histc(NE_i(:,2),linspace(min(NE_i(:)),max(NE_i(:)),100));
% c = histc(NE_i(:,3),linspace(min(NE_i(:)),max(NE_i(:)),100));
% figure,bar(edge,a,'facecolor','b')
% hold on, bar(edge,b,'facecolor','r')
% hold on, bar(edge,c,'facecolor','g')
% title('SNR Models: External Noise Inhibitory')
% legend('Stage1','Stage2','Stage3')

% Control hemisphere
edge = linspace(min(NE_C(:)),max(NE_C(:)),100);
a = histc(NE_C(:,1),linspace(min(NE_C(:)),max(NE_C(:)),100));
b = histc(NE_C(:,2),linspace(min(NE_C(:)),max(NE_C(:)),100));
c = histc(NE_C(:,3),linspace(min(NE_C(:)),max(NE_C(:)),100));
figure,bar(edge,a,'facecolor','b')
hold on, bar(edge,b,'facecolor','r')
hold on, bar(edge,c,'facecolor','g')
title('SNR Models CONTROL: External Noise Total')
legend('Stage1','Stage2','Stage3')

% edge = linspace(min(NE_C_e(:)),max(NE_C_e(:)),100);
% a = histc(NE_C_e(:,1),linspace(min(NE_C_e(:)),max(NE_C_e(:)),100));
% b = histc(NE_C_e(:,2),linspace(min(NE_C_e(:)),max(NE_C_e(:)),100));
% c = histc(NE_C_e(:,3),linspace(min(NE_C_e(:)),max(NE_C_e(:)),100));
% figure,bar(edge,a,'facecolor','b')
% hold on, bar(edge,b,'facecolor','r')
% hold on, bar(edge,c,'facecolor','g')
% title('SNR Models CONTROL: External Noise Excitatory')
% legend('Stage1','Stage2','Stage3')
% 
% edge = linspace(min(NE_C_i(:)),max(NE_C_i(:)),100);
% a = histc(NE_C_i(:,1),linspace(min(NE_C_i(:)),max(NE_C_i(:)),100));
% b = histc(NE_C_i(:,2),linspace(min(NE_C_i(:)),max(NE_C_i(:)),100));
% c = histc(NE_C_i(:,3),linspace(min(NE_C_i(:)),max(NE_C_i(:)),100));
% figure,bar(edge,a,'facecolor','b')
% hold on, bar(edge,b,'facecolor','r')
% hold on, bar(edge,c,'facecolor','g')
% title('SNR Models CONTROL: External Noise Inhibitory')
% legend('Stage1','Stage2','Stage3')
