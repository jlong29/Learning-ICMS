clear
cd('c:\program files\matlab\r2006b\work\Learning_ICMS_figures')
load('Acc_and_Spike_example.mat')

scrnsz = get(0,'screensize');
figure('Position',scrnsz)

% Lowpass filter Acc data at 50 Hz
[B,A] = butter(5,50/1000,'low');
Acc= filter(B,A,Acc);

plot(time_line,HF_Spikes(:,1)+10,'k'),hold on
% plot(time_line,HF_Spikes(:,2)+8,'k')
% plot(time_line,HF_Spikes(:,3)+6,'k')
plot(time_line,HF_Spikes(:,4)+8,'k')
plot(time_line,HF_Spikes(:,5)+6,'k')

plot(time_line2,Acc(:,1).*1.5+4,'b','linewidth',2)
plot(time_line2,Acc(:,3).*1.5+1.5,'b','linewidth',2)
plot(time_line2,Acc(:,2).*1.5-2,'b','linewidth',2)