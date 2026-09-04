function lim=datalim(data)
%datalim
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to calculate an appropriate limit for data, so that for better figure plot of imagesc function
if nargin==1
    s=0.01;
end
data(isnan(data))=0;
data(abs(data)==Inf)=0;
a=data(data~=0);
a1=unique(a);
n=length(a1);
n1=round(n*s);
lim=[a1(max(1,n1)),a1(max(1,n-n1))];
end

