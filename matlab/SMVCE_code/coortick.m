function [xtl,ytl,xtick,ytick]=coortick(coor,dxy)
% set the xy ticklabel as the lonlat based on the coor
% the returned xtl ytl is the xticklabel and yticklabel
if ~exist('dxy','var')
    dxy=[1,1];
end
dx=coor.post_lon*coor.width;
dy=-coor.post_lat*coor.nlines;
Xs=floor(coor.corner_lon);
while Xs<coor.corner_lon
    Xs=Xs+dxy(1);
end
XX=Xs:dxy(1):coor.corner_lon+dx;

Ys=floor(coor.corner_lat-dy);
while Ys<coor.corner_lat-dy
    Ys=Ys+dxy(2);
end
YY=Ys:dxy(2):coor.corner_lat;

xy1x=[XX',XX'];
xy1y=[YY',YY'];

xytickx=lonlat2sub(xy1x,coor);
xyticky=lonlat2sub(xy1y,coor);
xtick=xytickx(:,2);
ytick=xyticky(end:-1:1,1);
n=size(xy1y,1);
% xtl=cell(n,1);
ytl=cell(n,1);
for i=1:n
%     xtl{i}=num2str(round(xy1(i,1)*1000)/1000);
    ytl{i}=[num2str(round(xy1y(n-i+1,2)*1000)/1000),'{\circ}N'];
end

n=size(xy1x,1);
xtl=cell(n,1);
% xt=zeros(1,n);
% ytl=xtl;
for i=1:n
    xtl{i}=[num2str(round(xy1x(i,1)*1000)/1000),'{\circ}E'];
%     ytl{i}=num2str(round(xy1(i,2)*1000)/1000);
end

set(gca,'xticklabel',xtl,'xtick',xtick,'ytick',ytick,'yticklabel',ytl,'YTickLabelRotation',90);
end
function lonlat=sub2lonlat(subs,coor)
lonlat=[(subs(:,2)-1)*coor.post_lon+coor.corner_lon,...
    (subs(:,1)-1)*coor.post_lat+coor.corner_lat];
end
function [subs,subs1]=lonlat2sub(lonlat,coor)
%get the sub index of lonlat based on the coordinate of coor
subs=[round((lonlat(:,2)-coor.corner_lat)/coor.post_lat+1),...
    round((lonlat(:,1)-coor.corner_lon)/coor.post_lon+1)];
subs1=[((lonlat(:,2)-coor.corner_lat)/coor.post_lat+1),...
    ((lonlat(:,1)-coor.corner_lon)/coor.post_lon+1)];
end