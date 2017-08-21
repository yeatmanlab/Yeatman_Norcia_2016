basedir='~/git/Yeatman_Norcia_2016/data/'
subs={'TxtFaceN=11_20110906'};

subdir=fullfile(basedir,subs{1},'Exp_MATL_HCN_128_Avg');
% parameters for calculating implicit time
condt={'001' '002' '003' '004' '005' '006' '007'};
condf={'012' '013' '014' '015' '016' '017' '018'};
pcut=[.05 3]; weightFit=0; makeFigs=0; har=1;
% Calculate implicit time and amplitude for text conditions
[dTt degt degunwt radt hzt yVt yVerrt degt_se dTt_se] = mrC_CalcPhaseFreqSlope(subdir,condt,pcut,weightFit,makeFigs,har);
% Calculate implicit time and amplitude for face conditions
[dTf degf degunwf radf hzf yVf yVerrf degf_se dTf_se] = mrC_CalcPhaseFreqSlope(subdir,condf,pcut,weightFit,makeFigs,har);
% Make phase frequency plots for text
enums = [91 68]
c = vals2colormap(dTt,'parula',[140 260]);
figure;hold
for ii = enums
    errorbar(hzt(1:5),degunwt(1:5,ii),degt_se(1:5,ii),'o','color',c(ii,:),'markerfacecolor',c(ii,:));
    plot(hzt(1:5),polyval(polyfit(hzt(1:5)',degunwt(1:5,ii),1),hzt(1:5)),'-','color',c(ii,:))
end
axis tight
dTt(enums)
dTt_se(enums)
% Make phase frequency plots for face
enums = [82 120]
c = vals2colormap(dTf,'parula',[140 260]);
figure;hold
for ii = enums
    errorbar(hzf(1:5),degunwf(1:5,ii),degf_se(1:5,ii),'o','color',c(ii,:),'markerfacecolor',c(ii,:));
    plot(hzf(1:5),polyval(polyfit(hzf(1:5)',degunwf(1:5,ii),1),hzf(1:5)),'-','color',c(ii,:))
end
axis tight
dTf(enums)
dTf_se(enums)

% Find the electrodes with max signal across frequencies
mVt=mean(yVt);
mVf=mean(yVf);
[Vtsort, Vtind]=sort(mVt,2,'descend');
[Vfsort, Vfind]=sort(mVf,2,'descend');

% Plot amplitude as a function of frequency for the two electrodes with max
% signal
% First for the word electrode
figure;hold
nelec = 1;%[1:4];
errorbar(hzt,mean(yVt(:,Vtind(nelec)),2),mean(yVerrt(:,Vtind(nelec)),2),'-o',...
    'color', [.6 .6 .65],'markerfacecolor',[.6 .6 .65],'markersize',10,'linewidth',2);
errorbar(hzt,mean(yVf(:,Vtind(nelec)),2),mean(yVerrf(:,Vtind(nelec)),2),'-o',...
    'color', [0 0 0],'markerfacecolor',[0 0 0],'markersize',10,'linewidth',2);
set(gca,'xtick',hzt,'xticklabel',hzt)
axis([0 13 0 3.1])
% First for the face electrode
figure;hold
errorbar(hzt,mean(yVt(:,Vfind(nelec)),2),mean(yVerrt(:,Vfind(nelec)),2),'-o',...
    'color', [.6 .6 .65],'markerfacecolor',[.6 .6 .65],'markersize',10,'linewidth',2);
errorbar(hzt,mean(yVf(:,Vfind(nelec)),2),mean(yVerrf(:,Vfind(nelec)),2),'-o',...
    'color', [0 0 0],'markerfacecolor',[0 0 0],'markersize',10,'linewidth',2);
set(gca,'xtick',hzt,'xticklabel',hzt)
axis([0 13 0 3.1])

%% Repeat for even harmonics
subdir=fullfile(basedir,subs{1},'Exp_MATL_HCN_128_Avg');
% parameters for calculating implicit time
condt={'008' '009' '010' '011'}; 
condf={'019' '020' '021' '022'};
pcut=[.05 3]; weightFit=0; makeFigs=1; har=2;
% Calculate implicit time and amplitude for text conditions
[dTt degt degunwt radt hzt yVt yVerrt] = mrC_CalcPhaseFreqSlope(subdir,condt,pcut,weightFit,makeFigs,har);
% Calculate implicit time and amplitude for face conditions
[dTf degf degunwf radf hzf yVf yVerrf] = mrC_CalcPhaseFreqSlope(subdir,condf,pcut,weightFit,makeFigs,har);
% % Find the electrodes with max signal across frequencies
mVt=mean(yVt);
mVf=mean(yVf);
[Vtsort, Vtind]=sort(mVt,2,'descend');
[Vfsort, Vfind]=sort(mVf,2,'descend');

% Plot amplitude as a function of frequency for the two electrodes with max
% signal
% First for the word electrode
figure;hold
errorbar(hzt,mean(yVt(:,Vtind(1:4)),2),mean(yVerrt(:,Vtind(1:4)),2),'-o',...
    'color', [.6 .6 .65],'markerfacecolor',[.6 .6 .65],'markersize',10,'linewidth',2);
errorbar(hzt,mean(yVf(:,Vtind(1:4)),2),mean(yVerrf(:,Vtind(1:4)),2),'-o',...
    'color', [0 0 0],'markerfacecolor',[0 0 0],'markersize',10,'linewidth',2);
set(gca,'xtick',hzt,'xticklabel',hzt)
% First for the face electrode
figure;hold
errorbar(hzt,mean(yVt(:,Vfind(1:4)),2),mean(yVerrt(:,Vfind(1:4)),2),'-o',...
    'color', [.6 .6 .65],'markerfacecolor',[.6 .6 .65],'markersize',10,'linewidth',2);
errorbar(hzt,mean(yVf(:,Vfind(1:4)),2),mean(yVerrf(:,Vfind(1:4)),2),'-o',...
    'color', [0 0 0],'markerfacecolor',[0 0 0],'markersize',10,'linewidth',2);
set(gca,'xtick',hzt,'xticklabel',hzt)


