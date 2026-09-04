function lonlat=sub2lonlat(subs,coor)
lonlat=[(subs(:,2)-1)*coor.post_lon+coor.corner_lon,...
    (subs(:,1)-1)*coor.post_lat+coor.corner_lat];
end