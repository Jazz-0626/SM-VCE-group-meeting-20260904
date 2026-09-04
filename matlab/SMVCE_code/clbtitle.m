function []=clbtitle(clbh,str,fs)
%clbtitle
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%--This function is used to set a title of a colorbar, for example, the unit of the figure

% clbh=colorbar;
% str='m';
if ~exist('fs','var')
    fs=10;
end
h=get(clbh,'title');
set(h,'string',str,'fontsize',fs);

pos0=get(clbh,'position');
% x1=pos0(1)+pos0(3)/2;
% y1=pos0(2)+pos0(4)/2;

set(h,'units','normalized');
set(h,'position',[0.5,0.5,0]);
set(h,'horizontalalignment','center','verticalalignment','middle');

if pos0(3)<pos0(4)
    set(get(clbh,'title'),'rotation',90);
end
end
