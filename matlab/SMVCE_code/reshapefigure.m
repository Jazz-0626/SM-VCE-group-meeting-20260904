function reshapefigure(coor,axsi)

%the ratio between height and width
dis_row=110;
dis_col=110*cosd(coor.corner_lat+coor.post_lat*coor.nlines/2);
sca=-dis_row*coor.nlines*coor.post_lat/(dis_col*coor.width*coor.post_lon);

posfigure=get(gcf,'Position');
if ~exist('axsi','var')
posi=get(gca,'Position');
else
    axs0=gca;
    axes(axsi);
    hold on
    posi=get(gca,'Position');

    axes(axs0);
    hold on
end
scai=posi(4)*posfigure(4)/(posi(3)*posfigure(3));
posfigure(4)=posfigure(4)*sca/scai;
set(gcf,'Position',posfigure);

end