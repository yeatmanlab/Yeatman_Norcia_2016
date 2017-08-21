for ii=1:size(yV,1)
    figure;hold;
    title(sprintf('1F1 Amplitude %s HZ',num2str(hz(ii))),'fontsize',18,'fontname','times');
    bar(yV(ii,:),'facecolor','r');
    axis([0 128 0 max(yV(:))])
    errorbar(1:128,yV(ii,:),yVerr(ii,:), '.k','markersize',8);
    ylabel('Amplitude uV','fontsize',18,'fontname','times');
    xlabel('Electrode','fontsize',18,'fontname','times');
    set(gca,'fontsize',18,'fontname','times');
    figure;plotOnEgi(yV(ii,:));cbid=colorbar;colormap jet;caxis([0 max(yV(:))]);
    set(gca,'xtick',[],'ytick',[]);
    set(get(cbid,'ylabel'),'String','Amplitude uV');
    title(sprintf('1F1 Amplitude %s HZ',num2str(hz(ii))),'fontsize',14,'fontname','times');;
end