function [sita1,P_vce,x,f,iterator]=SMVCE_vce(L,B,P,enanglei)
%SMVCE_vce: determining the weight of observations L based on Variance Component Estimation
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
% Input:
%   L,B,P: each is a 1*N cell (N is the number of observation types), L1=B1*x, P1
%           P is not madatary.
% Output:
%   sita1:  the variance of unit weight, N*1
%   P_vce:  same as the input P, but determined by the VCE
%   x:      the solved unknown vector
%   f:      a flag indicate the VCE result, 0 is ok, 1 indicates appearing
%               negtive sita, 2 indicates exceeding the maximum iteration
%   iterator:the iteration number of VCE
%
%
vce_mode=0;
t1=2;
t2=5;
thr1=0.99;
thr2=0.0005;
thr3=0.000;
data_num0=length(L);
k=zeros(data_num0,1);
f=0;

for i=1:data_num0
    k(i)=length(L{i});
end
if nargin==2
    P=cell(1,data_num0);
    for i=1:data_num0
        P{i}=sparse(eye(k(i)));
    end
end
for i=1:data_num0
    if k(i)/max(k)<thr3
        L{i}=[];
        P{i}(:)=0;
    end
end

P0=P;
[sita,v,x]=getsita(B,P,L,enanglei);
k(k/max(k)<thr3)=0;
neq0=(find(k~=0))';

pscales=diag(P{neq0(1)});pscales(pscales==0)=[];
pscale=mean(pscales);
iterator=0;
while (min(sita)/max(sita)<thr1&&max(sita)-min(sita)>thr2)||((sum(sita<0))>0)
    iterator=iterator+1;
    if (sum(sita<=0))>0 %||max(sita)/min(sita)>10^10%%%%%%%%%%%%%%%%%%%%%%%%%% ################
        if (sum(sita>0))<3
            sita1=zeros(length(L),1);
            P_vce=P0;
            f=1;
            return;
        else
            iterator=1;
            P=P0;
            idxi=neq0(sita<=0);
            for idxii=idxi
                P{idxii}(:)=0;
                L{idxii}=[];
            end
            neq0(sita<=0)=[];
            [sita,v]=getsita(B,P,L,enanglei);
        end
    end
    if iterator>5
        P_vce=P;
        sita1=zeros(length(L),1);
        sita1(neq0)=sita;
        f=2;
        return;
    end
    
    [~,P]=getP(sita,v,P,neq0,t1,t2,vce_mode);
    [sita,v,x]=getsita(B,P,L,enanglei);
    
    pscales=diag(P{neq0(1)});pscales(pscales==0)=[];
    pscale=mean(pscales);
end

P_vce=P;
sita1=zeros(length(L),1);
sita1(neq0)=sita;
end

%% other functions
function [sita,v,x]=getsita(B,P,L,enanglei)
data_num0=length(L);
k=zeros(data_num0,1);
unknum=size(B{1},2);
N=zeros(unknum,unknum);
U=zeros(unknum,1);
for i=1:data_num0
    k(i)=length(L{i});
end
%     k(k<sqrt(max(k)))=0;
% k(k<10)=0;
neq0=(find(k~=0))';
data_num=length(neq0);
k=k(k~=0);
Ni=cell(data_num,1);
v=cell(data_num,1);
ii=0;
for i=neq0
    ii=ii+1;
    P{i}(isnan(P{i}))=0;
    Ni{ii}=B{i}'*P{i}*B{i};
    N=N+Ni{ii};
    U=U+B{i}'*P{i}*L{i};
end
N(abs(N)==inf)=nan;
if sum(isnan(N(:)))>0
    sita=zeros(1,length(neq0));
    v=0;
    x=zeros(unknum,1);
    return;
end


        tanv=tan(enanglei);
        Bv=[tanv,-1,0]; 
        tanvs=3;
        ftanvs=0;
        if abs(tanv)>tanvs
            ftanvs=1;
           Bv=[tanv/abs(tanv)*tanvs,-1,0]; 
        end        
        Nv=Bv'*Bv;
        N1=N;
        N1(1:3,1:3)=N1(1:3,1:3)+Nv;
        
NN=pinv(N1);
x=NN*U;

if sign(x(2)*sin(enanglei))<0%&&ftanvs==1
tanv=-tan(enanglei);
        Bv=[tanv,-1,0]; 
        tanvs=3;
        if abs(tanv)>tanvs
           Bv=[tanv/abs(tanv)*tanvs,-1,0]; 
        end        
        Nv=Bv'*Bv;
        N1=N;
        N1(1:3,1:3)=N1(1:3,1:3)+Nv;
        
NN=pinv(N1);
x=NN*U;
end
W=zeros(data_num,1);
ii=0;
for i=neq0
    ii=ii+1;
    v{ii}=B{i}*x-L{i};
    W(ii)=v{ii}'*P{i}*v{ii};
end

S=zeros(data_num,data_num);
for i=1:data_num
    for j=1:data_num
        if i==j
            S(i,j)=k(i)-2*trace(NN*Ni{i})+trace(NN*Ni{i}*NN*Ni{i});
            continue;
        end
        S(i,j)=trace(NN*Ni{i}*(NN*Ni{j}));
    end
end
sita=(pinv(S)*W);
end
function P0=getP0(v,P,neq0,t2)
ii=0;
P0=P;
for i=neq0
    ii=ii+1;
    v1=v{ii};
    p1=diag(P{i});
    v1(p1==0)=nan;
    avg=nanmedian(v1);
    v2=[abs(v1-avg),[1:length(v1)]'];
    v3=sortrows(v2,1);
    sd=nanstd(v1(v3(1:round(0.68*length(v1)),2)))*2;
    b=abs(v1-avg)/sd;
    a=find(b>t2);
    p1(b>t2)=0;
    P0{i}=diag(p1);
end
end
function [P0,P]=getP(sita,v,P,neq0,t1,t2,vce_mode)
ii=0;
P0=P;
for i=neq0
    ii=ii+1;
    P{i}=sita(1)/sita(ii)*P0{i};
    if vce_mode==1%&&iterator>0
        %2019/08/26
        P{i}=diag(robust_vce(diag(P{i}),diag(P0{i}),v{ii},sita(ii),t1,t2));
%         P{i}=robust_vce1(P{i},P0{i},v{ii},sita(ii));
    end
end
end

function y=robust_vce(P_vec,P_vec0,v,sita0,t1,t2)
Pi_tem=P_vec;
P_vec0(P_vec0==0)=[];
b=v;
d_tem=sita0/mean(P_vec0);
d_tem=(median(abs(v))/0.6745)^2;
b1=abs(b/sqrt(d_tem));
a=find(abs(b1)>=t2);
Pi_tem(abs(b1)>=t2)=0;
t1t21=find(((b1<t2).*(b1>=t1))'==1);
Pi_tem(t1t21)=t1*Pi_tem(t1t21)./b1(t1t21).*(((t2-abs(b1(t1t21)))/(t2-t1)).^2);
y=Pi_tem;
end

function y=robust_vce1(P_vec,P_vec0,v,sita)
ps=diag(P_vec0);ps(ps==0)=[];
sd=sita/mean(ps);
R=v/(2.385*sd);
W=1./(1+R.^2);
W=W/mean(W);
y=P_vec.*diag(W);
end
