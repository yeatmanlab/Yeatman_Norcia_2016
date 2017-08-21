datadir=pwd
cond={'001' '002' '003' '004' '005'};
pcut=[.05 4];
weightFit=0;
makeFigs=0;
har=[];
[dT deg degunw rad hz yV yVerr sedeg sdT] = mrC_CalcPhaseFreqSlope(datadir,cond,pcut,weightFit,makeFigs,har);
figure;
errorbar(hz,degunw(:,59),sedeg(:,59),'ko');
lsline;
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
axis square
axis([0 7 0 550])
title(sprintf('Text, Electrode 59, Implicit time = %.2fms +/- %.2fms',dT(59),sdT(59)))
%Faces
cond={'012' '013' '014' '015' '016'}
[dT deg degunw rad hz yV yVerr sedeg sdT] = mrC_CalcPhaseFreqSlope(datadir,cond,pcut,weightFit,makeFigs,har);
figure;
errorbar(hz,degunw(:,91),sedeg(:,91),'ko');
lsline;
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')
title(sprintf('Faces, Electrode 91, Implicit time = %.2fms +/- %.2fms',dT(91),sdT(91))) 
axis square
axis([0 7 0 500])