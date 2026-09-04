function strs=splitstr(str)
%splitstr
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to split a string into several parts
seperators=[',',':',' ','	'];
flag=ismember(str,seperators);
inds=find(flag==0);
di=inds(2:end)-inds(1:end-1);
inds2=find(di>1);
inds2=[0,inds2,length(inds)];
strs=cell(1,length(inds2)-1);
if ~isempty(inds2)
    for i=1:length(strs)
        indsii=inds(inds2(i)+1):inds(inds2(i+1));
        strs{i}=str(indsii);
    end
end
end
