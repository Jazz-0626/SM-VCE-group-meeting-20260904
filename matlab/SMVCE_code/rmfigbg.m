function []=rmfigbg(pfig)
%rmfigbg
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to remove the blank margin of a grid picture
fig=imread(pfig);
a=sum(fig==255,3);
b=sum(fig==0,3);
b(a==3)=3;
b1=find(sum(b~=3,1)~=0);
b2=find(sum(b~=3,2)~=0);
fig1=fig(min(b2):max(b2),min(b1):max(b1),:);
imwrite(fig1,pfig);
end
