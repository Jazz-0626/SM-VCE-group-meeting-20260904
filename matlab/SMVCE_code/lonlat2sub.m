function subs=lonlat2sub(lonlat,coor)
%lonlat2sub
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to get the sub index of lonlat based on the coordinate of coor
subs=[round((lonlat(:,2)-coor.corner_lat)/coor.post_lat+1),...
    round((lonlat(:,1)-coor.corner_lon)/coor.post_lon+1)];
end
