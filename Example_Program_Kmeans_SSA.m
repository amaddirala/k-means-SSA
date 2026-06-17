clc;
clear all;
close all;
addpath('FastICA_25\')
load('ERP_EEG1.mat');
load('ERP_EOG1.mat');
fs=256;
x = s+a;time = [0:length(x)-1]/fs;
M = 128;clusters = 4;Th = 1.4; Tssa = 0.01;nstd=0.1;ne=100;levels = 4;order = 8;threshold = 500;

figure;
[eemd_ica_eog]=My_EEMD_ICA(x,nstd,ne);%% you have to enter the index of eye-blink component EEMD-ICA
[ssa_ica_eog]=My_SSA_ICA(x,M,levels,fs,order);%% you have to enter the index of eye-blink component SSA-ICA
[shenghuan_eog] = Shenghuan_technique(x,threshold,fs);
[fbse_ewt_eog]=My_FBSE_EWT(x,fs,fs,fs/4);
[proposed_eog]=My_kmeans_SSA_Method(x,M,clusters,Th,Tssa);
close all;
figure;
subplot(6,2,1);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,a,'Color','k');hold on;box off

subplot(6,2,2);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,s,'Color',[0.3 0.75 0.93]);hold on;box off

subplot(6,2,3);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,eemd_ica_eog,'Color','b');hold on;box off;title(['EEMD-ICA ($\hat{a}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,4);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,x - eemd_ica_eog,'Color','b');hold on;box off;title(['EEMD-ICA ($\hat{s}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,5);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,ssa_ica_eog,'Color','b');hold on;box off;title(['SSA-ICA ($\hat{a}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,6);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,x - ssa_ica_eog,'Color','b');hold on;box off;title(['EEMD-ICA ($\hat{s}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,7);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,shenghuan_eog,'Color','b');hold on;box off;title(['Method[45] ($\hat{a}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,8);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,x - shenghuan_eog,'Color','b');hold on;box off;title(['Method[45] ($\hat{s}$)'],'FontWeight','normal','Interpreter','latex')


subplot(6,2,9);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,fbse_ewt_eog,'Color','b');hold on;box off;title(['FBSE-EWT ($\hat{a}$)'],'FontWeight','normal','Interpreter','latex')


subplot(6,2,10);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,x - fbse_ewt_eog,'Color','b');hold on;box off;title(['FBSE-EWT ($\hat{s}$)'],'FontWeight','normal','Interpreter','latex')

subplot(6,2,11);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,proposed_eog,'Color','b');hold on;box off;title(['Proposed ($\hat{a}$)'],'FontWeight','normal','Interpreter','latex')
xlabel('Time (s)')
subplot(6,2,12);plot(time,x,'Color',[0.93 0.69 0.13]);hold on;
plot(time,x - proposed_eog,'Color','b');hold on;box off;title(['Proposed ($\hat{s}$)'],'FontWeight','normal','Interpreter','latex')
xlabel('Time (s)')
set(gcf, 'Position',[100 100 800 600]);