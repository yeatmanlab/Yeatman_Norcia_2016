function [dT deg degunw rad hz yV yVerr sedeg sdT] = mrC_CalcPhaseFreqSlope(datadir,cond,pcut,weightFit,makeFigs,har)
% Copmute time delay based on the slope of change over frequency.
%
% Inputs   datadir = The directory with the exported power diva session
%          cond    = A cell aray with condition numbers that you would like
%                    to use to compute phase/frequency slope
%          pcut    = The minimum pvalue for which to use a phase
%                    estimate in the linear fit and also how many
%                    significant frequencies there must be to fit a line
%          weightFit Weight the fitting of the slope by gaussian distance.
%                    The input will be the maximum amount of weighting
%          har     = Which harmonic to use 1st 2nd 3rd
%          makeFigs= Make figures? 0 or 1
% Outputs  deg     = Phase in degrees
%          degunw  = Phase in degrees unwrapped
%          rad     = Phase in radians
%          dT      = Time delay in miliseconds
%          hz      = The hz of the stimulus
%
%
% Example (for text face study):
% %For Text
% cond={'001' '002' '003' '004' '005'};
% datadir= '/Users/jyeatman/Desktop/Exp_MATL_HCN_128_Avg';
% [deg degunw rad dT dTCI] = mrC_CalcPhaseFreqSlope(datadir,cond);
%
% %For faces
% cond={'012' '013' '014' '015' '016'}

if ~exist('pcut','var') || isempty(pcut)
    pcut = [0.05 4];
end
if ~exist('makeFigs','var') || isempty(makeFigs)
    makeFigs=0;
end
if ~exist('weightFit','var') || isempty(weightFit)
    weightFit=0;
end
if ~exist('har','var') || isempty(har)
    har=1; %default to first harmonic
end
%% loop over the number of conditions and collect data
for ii = 1 : length(cond)
    % Load condition file
    condFile=fullfile(datadir,['Axx_c' cond{ii} '.mat']);
    load(condFile);
    % The row in the data with the 1F1
    F1(ii)=har*i1F1+1;
    % Frequency in hz for this condition
    hz(ii)=i1F1./2;
    % Amplitude
    yV(ii,:) = Amp(F1(ii),:);
    yVerr(ii,:)=SpecStdErr(F1(ii),:);
    % Calculate phase in radians for the 1F1
    rad(ii,:) = atan2(Sin(F1(ii),:),Cos(F1(ii),:));
    % Get the pvalue for the amplitude/phase estimates and calculate a
    % zscore
    pval(ii,:) = SpecPValue(F1(ii),:);
    z(ii,:) = norminv(1-pval(ii,:));
    % Calculate the standard error on the phase estimate
    serad(ii,:)=atan2(SpecStdErr(F1(ii),:),Amp(F1(ii),:));
end
% Determine maximum zscore for weighting fits
z(z > weightFit)=weightFit;
% The interval of atan2 is from -pi to pi. Convert to 0 to 2pi, meaning
% that -pi/2 is converted to 3pi/2
neg = rad<0;
rad(neg) = 2*pi + rad(neg);
% Convert radians to degrees
deg = rad.*(180/pi);
sedeg=serad.*(180/pi);
%% Phase unwrap
degunw=deg;
for ii = 2 : size(deg,1)
    % check if any conditions phase estimate is less than the previous
    % condition
    while sum(degunw(ii,:)<degunw(ii-1,:))>0
        % add 360 to any phase that is less than the previous one
        unw=(degunw(ii,:)<degunw(ii-1,:)).*360;
        degunw(ii,:)=degunw(ii,:)+unw;
        clear unw
    end
end
%% Fit a line to frequency conditions
plotnum=0;
b=zeros(length(degunw),2);
sdFit=zeros(length(degunw),2);
for ii = 1 : length(degunw)
    % only use significant phase estimates
    use = pval(:,ii) < pcut(1);
    if sum(use)>=pcut(2)
        phase = degunw(use,ii);
        sephase = sedeg(use,ii);
        freq = hz(use)';
        if weightFit>0
            % Weight the fit by the zscore.  Multiplying by 100 means that
            % rounding should not cause problems
            phaseW=[];
            freqW=[];
            for jj=1:length(phase)
                phaseW(length(phaseW)+1:length(phaseW)+round(100*z(jj,ii)))=repmat(phase(jj),round(100*z(jj,ii)),1);
                freqW(length(freqW)+1:length(freqW)+round(100*z(jj,ii)))=repmat(freq(jj),round(100*z(jj,ii)),1);
            end
            b(ii,:) = polyfit(freqW ,phaseW,1);
            %bootsrap fits to get an error measure
            for k=1:500
                num=randn(sum(use),1).*sephase;
                phasek=phase+num;
                bootstrp(k,:)=polyfit(freq,phasek,1);
                clear phasek num
            end
            sdFit(ii,:)=std(bootstrp);
        else
            b(ii,:) = polyfit(freq ,phase,1);
            for k=1:500
                num=randn(sum(use),1).*sephase;
                phasek=phase+num;
                bootstrp(k,:)=polyfit(freq,phasek,1);
                clear phasek num
            end
            sdFit(ii,:)=std(bootstrp);
        end
        
    end
    %% Plot the fits if desired
    if  makeFigs==1
        if sum(ii==[1 26 51 76 101 126])==1
            figure;
            plotnum=0;
        end
        plotnum=plotnum+1;
        subplot(5,5,plotnum);hold;
        plot(hz(use)',degunw(use,ii),'b.');
        plot(hz(~use)',degunw(~use,ii),'r.');
        plot([0 hz(end)]',polyval(b(ii,:),[0 hz(end)])','-k');
        title(sprintf('Electrode %s',num2str(ii)));
        set(gca,'xtick',[2:2:hz(end)],'ytick',[360:360:round(max(degunw(pval<pcut(1))))],'fontsize',8,'fontname','times');
        axis([0 hz(end) 0 round(max(degunw(pval<pcut(1))))],'square');
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % %
    % %This is code to fit the function with polyfit instead of regress. It may
    % %be more flexible if we ever want to fit curves.
    %     s(ii,:) = polyfit(hz ,degunw(:,ii)',1);
    %     yfit = polyval(s(ii,:),hz);
    %     yresid = degunw(:,ii)' - yfit;
    %     SSresid = sum(yresid.^2);
    %     SStotal = length(degunw(:,ii)'-1) * var(degunw(:,ii)');
    %     R2(ii) = 1- SSresid/SStotal;
    %     clear yfit yresid SSresid SStotal
    % % % % % % % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % %
end

% Calculate time lag from slope
dT=b(:,1)./0.36;
sdT=sdFit(:,1)./0.36;
% Make a figure of time lage for each channnel
figure;
bar(dT, 'facecolor',[1 .5 0]);
hold;
m=median(dT(dT>0));
plot([1 129],[m m],'-k');
axis([0 129 0 400]);
set(gca,'xtick',1:3:128);
xlabel('Electrode');ylabel('Time Delay');
dT(dT==0)=nan;
sdT(sdT==0)=nan;
% To add error bars
x=1:128;
errorbar(x,dT,sdT, '.k','markersize',4)
hold;
plotDT=dT;
%plotDT(dT>250)=240;
%plotDT(dT<150)=160;
figure;plotOnEgi(plotDT);cbid=colorbar;colormap jet;
title('Implicit time');
%set(gca,'xtick',[],'ytick',[]);
%set(get(cbid,'ylabel'),'String','Implicit Time')

%% Make amplitude plots
figure;
ymax = max(yV(:));
for ii=1:length(hz)
    subplot(1,length(hz),ii);hold;
    bar(yV(ii,:), 'facecolor',[1 .5 0]);
    plot([1:128;1:128],[yV(ii,:)-yVerr(ii,:);yV(ii,:)+yVerr(ii,:)],'-k');
    ylabel('Amplitude yV');
    title(sprintf('%d Hz',hz(ii)))
    xlim([0 129])
    ylim([0 ymax]);
end

% Set colorbar to go from 20 to 80th percentiles
prc = prctile(yV(:),[10 90]);
for ii=1:length(hz)
    figure;
    plotOnEgi(yV(ii,:));colorbar;colormap jet;
    caxis(prc);
    title(sprintf('Amplitude (yV) %d Hz',hz(ii)));
end
return
%
% % To check what the HZ is for each condition in an exported session.
% d=dir('Axx*.mat')
% for ii=1:length(d)
%     load(d(ii).name);
%     % The row in the data with the 1F1
%     F1(ii)=i1F1+1;
%     hz(ii)=i1F1./2;
%     oddHarmonic(ii,:)=Amp(F1(ii),:)
% end