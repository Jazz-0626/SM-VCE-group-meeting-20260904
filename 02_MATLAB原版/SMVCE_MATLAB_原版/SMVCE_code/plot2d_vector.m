function xyend=plot2d_vector(data,scale,arrorcolor)
% By Jihong Liu @ KAUST
% plot the 2d_vector
% Input:
%   data:Nps*4, lon,lat,east,north
%   scale: to indicate the length of the arrow line, default:1
%   arrorcolor: color of arrows

%%%%%%%%%%% Arrow size %%%%%%%%%%%%%%%%%%
K_arrow=0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlim=get(gca,'xlim');
ylim=get(gca,'ylim');

ydir=get(gca,'ydir');
if strcmp(ydir,'reverse')
    data(:,4)=-data(:,4);
end
xdir=get(gca,'xdir');
if strcmp(xdir,'reverse')
    data(:,3)=-data(:,3);
end
if ~exist('scale','var')||scale==0
    scale=1;
end
if ~exist('lw','var')||lw==0
    lw=1;
end
if ~exist('fs','var')||fs==0
    fs=15;
end
if ~exist('arrorcolor','var')
    arrorcolor=[0,0.45,0.74];
end

Nps=size(data,1);
theta= pi / 12;% 箭头角度
A1 = [cos(theta), -sin(theta);sin(theta), cos(theta)]; % 旋转矩阵
theta = -theta;
A2 = [cos(theta), -sin(theta);sin(theta), cos(theta)];% 旋转矩阵
xyend=[];
for i=1:Nps
    xy0=data(i,1:2);
    dxy1=data(i,3:4)*scale;
    xy1=xy0+dxy1;
    xy=[xy0;xy1];
    hold on;
    arrow=-data(i,3:4)*K_arrow/sqrt(sum(data(i,3:4).^2))*scale/5;

    dxy1l=sqrt(sum(dxy1.^2));
    arrowl=sqrt(sum(arrow.^2));
    if arrowl>dxy1l
        arrow=arrow*(dxy1l/arrowl);
    end
    %     arrowSize=sqrt(sum(arrow.^2));
    %     arrow(arrow>=0)=arrowSize;
    %     arrow(arrow<0)=-arrowSize;
    arrow_1= A1 * arrow';
    arrow_2= A2 * arrow';
    arrow_1=  arrow_1 + xy1'; % 箭头的边的x坐标
    arrow_2=  arrow_2 + xy1'; % 箭头的变的y坐标
    % 三角箭头(填充)
    triangle_x= [arrow_1(1),xy1(1),arrow_2(1),arrow_1(1)];
    triangle_y= [arrow_1(2),xy1(2),arrow_2(2),arrow_1(2)];

    xy=[xy0;(arrow_1+arrow_2)'/2];
    plot(xy(:,1),xy(:,2),'-','linewidth',lw*2,'color',arrorcolor);
%     hi=fill(triangle_x,triangle_y,arrorcolor);
    hi=fill(triangle_x,triangle_y,arrorcolor);
    hi.LineWidth=0.5;
    xyend=[xyend;xy1];
end
if ~exist('data0','var')||length(data0)~=4
    xkongbai=0.25;
    ykongbai=0.1;
    if (exist('xlim','var')&&xlim(1)==0&&xlim(2)==1)&&(exist('ylim','var')&&ylim(1)==0&&ylim(2)==1)
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
    end
    medlen=ceil(nanmedian(sqrt(sum(data(:,3:4).^2,2))));
    data0=[xlim(2)-(xlim(2)-xlim(1))*xkongbai,ylim(1)+(ylim(2)-ylim(1))*ykongbai,medlen,0];
    if strcmp(ydir,'reverse')
        data0(2)=ylim(2)-(ylim(2)-ylim(1))*ykongbai;
    end
    if strcmp(xdir,'reverse')
        data0(1)=xlim(1)+(xlim(2)-xlim(1))*xkongbai;
    end
    if ~exist('str0','var')||~ischar(str0)
        str0=[num2str(medlen),'mm/yr'];
    else
        str0=[num2str(medlen),str0];
    end
end
xy0=data0(1:2);
xy1=xy0+data0(3:4)*scale;
xy=[xy0;xy1];
% text(xy1(1),xy1(2),str0,'edgecolor','white','fontsize',fs,'verticalalignment','middle','fontweight','bold','backgroundcolor','white');

arrow=-data0(3:4)*scale*K_arrow*2;
%     arrowSize=sqrt(sum(arrow.^2));
%     arrow(arrow>=0)=arrowSize;
%     arrow(arrow<0)=-arrowSize;
arrow_1= A1 * arrow';
arrow_2= A2 * arrow';
arrow_1=  arrow_1 + xy1'; % 箭头的边的x坐标
arrow_2=  arrow_2 + xy1'; % 箭头的变的y坐标
% 三角箭头(填充)
triangle_x= [xy0(1),xy0(1),xy1(1),xy1(1),xy0(1)];
triangle_y= [arrow_1(2),arrow_1(2),arrow_1(2),arrow_2(2),arrow_2(2)];
% fill(triangle_x,triangle_y,'white','edgecolor','white');

hold on;
% plot(xy(:,1),xy(:,2),'-','linewidth',lw*2,'color',arrorcolor);
triangle_x= [arrow_1(1),xy1(1),arrow_2(1),arrow_1(1)];
triangle_y= [arrow_1(2),xy1(2),arrow_2(2),arrow_1(2)];
% fill(triangle_x,triangle_y,arrorcolor);
end