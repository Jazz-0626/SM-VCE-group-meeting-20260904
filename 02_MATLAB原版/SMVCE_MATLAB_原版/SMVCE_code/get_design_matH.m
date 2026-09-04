function B=get_design_mat(Bgeo,de,dn,du,dims)
%get_design_mat
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to estalish the design matrix based on SM
%
%Inputs:
%  Bgeo: the geometry relationship between InSAR observations and 3-D
%        deformations, with size of m*n*(data_num*3), m,n is the size of
%        InSAR observations, data_num is the types of InSAR observations
%  de,dn,du: the coordinate increments with size of m*n with respect to the interested point
%Output:
%  B   : the design matrix based on SM
%
% Last modified: Sept. 20th, 2019

[m,n,data_num]=size(Bgeo);
data_num=data_num/3;
if mod(data_num,1)~=0
    fprintf('input Bgeo is incorrect...')
    return;
end

if nargin==4
    %considering 3-D SM
    B=zeros(m*n*data_num,12);
    for i=1:data_num
        ai=Bgeo(:,:,1+3*(i-1));
        bi=Bgeo(:,:,2+3*(i-1));
        ci=Bgeo(:,:,3+3*(i-1));
        a=ai(:);
        b=bi(:);
        c=ci(:);
        ae=ai(:).*de(:);
        an=ai(:).*dn(:);
        au=ai(:).*du(:);
        be=bi(:).*de(:);
        bn=bi(:).*dn(:);
        bu=bi(:).*du(:);
        ce=ci(:).*de(:);
        cn=ci(:).*dn(:);
        cu=ci(:).*du(:);
        B(1+m*n*(i-1):m*n*i,:)=[a,b,c,ae,an,au,be,bn,bu,ce,cn,cu];
    end
    ind=[];
    for i=1:3
        if ~ismember(i,dims)
            switch i
                case 1
                    ind=[ind,1,4,5,6];
                case 2
                    ind=[ind,2,7,8,9];
                case 3
                    ind=[ind,3,10,11,12];
            end
        end
    end
    B(:,ind)=[];
else
    %considering 2-D SM
    B=zeros(m*n*data_num,9);
    for i=1:data_num
        ai=Bgeo(:,:,1+3*(i-1));
        bi=Bgeo(:,:,2+3*(i-1));
        ci=Bgeo(:,:,3+3*(i-1));
        a=ai(:);
        b=bi(:);
        c=ci(:);
        ae=ai(:).*de(:);
        an=ai(:).*dn(:);
        be=bi(:).*de(:);
        bn=bi(:).*dn(:);
        ce=ci(:).*de(:);
        cn=ci(:).*dn(:);
        B(1+m*n*(i-1):m*n*i,:)=[a,b,c,ae,an,be,bn,ce,cn];
    end
    ind=[];
    for i=1:3
        if ~ismember(i,dims)
            switch i
                case 1
                    ind=[ind,1,4,5];
                case 2
                    ind=[ind,2,6,7];
                case 3
                    ind=[ind,3,8,9];
            end
        end
    end
    B(:,ind)=[];
end
end
