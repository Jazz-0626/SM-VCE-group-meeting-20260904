function lim=coorlim(coor)
%datalim
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%--This function is used to calculate the limit coor data
lim=[coor.corner_lon;
    coor.corner_lon+coor.post_lon*(coor.width-1);
    coor.corner_lat+coor.post_lat*(coor.nlines-1);
    coor.corner_lat];
end
