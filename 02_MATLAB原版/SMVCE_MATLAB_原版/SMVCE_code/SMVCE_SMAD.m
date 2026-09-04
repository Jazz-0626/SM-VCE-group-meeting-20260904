function [mask,smp3,defo0]=SMVCE_SMAD(smp,r0,c0,nom,iffault)
%SMAD
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to select points under the stationary assumption
%      based on the Strain Model dependent Adaptive Neighbourhood (SMAD)
%
%Input:
%  smp : a sample of data with size of m*n*data_num
%  r0,c0: the central point index
%  nom: L1/L2-norm for calculation
%  temps: the directional templates
%Output:
%  mask: a matrix with same size of smp, 1 / 0 indicates the cooresponding
%            pixels if or not included in the stationary assumption
%  smp3: smp3=smp.*mask;
%  defo0:the filtered deformation for the central point
%

[m,n,data_num]=size(smp);
remainratio=0.05;
defo0=zeros(1,data_num);
temps=gettemps([m,n],32);
if nargin<2
    r0=round(size(smp,1)/2);%the row index of the intersted point
end
if nargin<3
    c0=round(size(smp,2)/2);%the column index of the intersted point
end
if nargin<4
    nom=2;%solving mode, 1 L1-norm/2 L2-norm
end
if nargin<5
    iffault=0;%if there is fault trace data, there is no need for "getminsum"
end
if nargin<6
    sigscale1=3;%the threthod of residual/sigma in the first estimation
    sigscale2=6;%the threthod of residual/sigma in the second estimation
    Lratio1=0.65;%the threthod of residual/sigma in the first estimation
    Lratio2=1;%the threthod of residual/sigma in the second estimation
    flag_scd=1;% indicate if conduct the second estimation
    flag_show=0;% indicate if show the result
end
smp0=smp;
smp0(isnan(smp0))=0;
if sum(smp0(:)==0)==length(smp(:))
    mask=zeros(size(smp));
    smp3=zeros(size(smp));
    smp0=0;
    return;
end
%% Solving
%x->col;y->row
% ndir=size(temps,3);
mask=zeros(size(smp));
smp2=zeros(size(smp));
smp3=zeros(size(smp));
%inintal removing a half of pixels for estimating deformation gradients
%L=[1,dx,dy]*[L0,xg,yg]';
if iffault==0
    smp1=getminsum(smp,temps);
else
    smp1=smp;
end
% smp1=smp;
[x,y]=meshgrid(1:n,1:m);
x=x-c0;
y=y-r0;
figk=0;
dlsk=[];
for datai=1:data_num
    smpi=smp(:,:,datai);
    Length_smpi=length(smpi(smpi~=0));
    smp1i=smp1(:,:,datai);
    L0=smpi(:);
    L01=smp1i(:);
    B0=[ones(length(L01),1),x(:),y(:)];
    Lneq0=find(L01~=0);
    L02=L01(Lneq0);
    B02=B0(Lneq0,:);
    L02_med=[abs(L02-median(L02)),(1:length(L02))'];
    L02_med_sort=sortrows(L02_med,1);
    Lneq2=L02_med_sort(1:round(length(L02)*Lratio1),2);
    if length(Lneq2)<remainratio*Lratio1*Length_smpi
        continue;
    end
    L=L02(Lneq2);
    B=B02(Lneq2,:);
%     unk=(B'*B)\B'*L;
    unk=pinv(B'*B)*B'*L;
    rs1=B*unk-L;
    sig=sqrt(mean(rs1.^2));
    %     rs2=sort(rs1.^2);
    %the local variance of this data in this region,
    %use the 65% data to cal.
    %     sig=sqrt(mean(rs2(1:round(length(rs2)*0.65))));
    %     rs3=wls(rs2,ones(length(rs2),1),nom);
    %     sig=sqrt(rs3);
    %     if sig<0.001
    %         sig=0.01;
    %     end
    % calculate the residuals in the whole region
    rs0=B0*unk-smpi(:);
    maski=zeros(m,n);
    if sum(abs(rs0/(sig+eps))<=sigscale1)>remainratio*Length_smpi
        maski(abs(rs0/(sig+eps))<=sigscale1)=1;
    end
    
    smp2(:,:,datai)=smpi.*maski;
    %the second refinement
    if flag_scd==1&&sum(maski(:))~=length(maski(:))
        L02=L0(maski(:)==1);
        B02=B0(maski(:)==1,:);
        L02_med=[abs(L02-median(L02)),(1:length(L02))'];
        L02_med_sort=sortrows(L02_med,1);
        Lneq2=L02_med_sort(1:round(length(L02)*Lratio2),2);
        if length(Lneq2)<remainratio*Lratio2*Length_smpi
            defo0(datai)=unk(1);
            mask(:,:,datai)=maski;
            continue;
        end
        L=L02(Lneq2);
        B=B02(Lneq2,:);
%         unk=(B'*B)\B'*L;
        unk=pinv(B'*B)*B'*L;
        rs1=B*unk-L;
        sig=sqrt(mean(rs1.^2));
        %     rs3=wls(rs2,ones(length(rs2),1),nom);
        %     sig=sqrt(rs3);
        %         if sig<0.001
        %             sig=0.01;
        %         end
        % calculate the residuals in the whole region
        rs0=B0*unk-smpi(:);
        if sum(abs(rs0/(sig+eps))<=sigscale2)>remainratio*Length_smpi
            maski(abs(rs0/(sig+eps))<=sigscale2)=1;
        end
        smp3(:,:,datai)=smpi.*maski;
    else
        smp3(:,:,datai)=smp2(:,:,datai);
    end
    defo0(datai)=unk(1);
    mask(:,:,datai)=maski;
    if flag_show==1
        out(:,:,figk+1)=smp(:,:,datai);
        out(:,:,figk+2)=smp1(:,:,datai);
        out(:,:,figk+3)=reshape(rs0/sig,m,n);
        out(:,:,figk+4)=smp2(:,:,datai);
        out(:,:,figk+5)=smp3(:,:,datai);
        
        dl=[min(smpi(:)),max(smpi(:))];
        if sum(~isnan(smpi(:)))==0
            dl=[0,0];
        end
        if dl(1)==dl(2)
            dl=dl+[-0.5,0.5];
        end
        dls=[dl;dl;[0,0];dl;dl];
        dlsk=[dlsk;dls];
        figk=figk+5;
    end
end
% smp0=unk(1);
if flag_show==1
    if ~exist('getfig.m','file')
        fprintf('Note: not exist the function of getfig.m!\n');
    end
    xp=[c0-0.5,c0-0.5,c0+0.5,c0+0.5,c0-0.5];
    yp=[r0-0.5,r0+0.5,r0+0.5,r0-0.5,r0-0.5];
    out(out==0)=nan;
    out(abs(out)>10^10)=nan;
    axs=getfig(out,dlsk,1,0,[data_num,5]);
    for i=1:length(axs)
        axes(axs(i));hold on;
        plot(xp,yp,'-m','linewidth',2);
    end
end
end
%% other functions
function smp1=getminsum(smp,temps)
data_num=size(smp,3);
ndir=size(temps,3);
smp(smp==0)=nan;
% temps=gettemps([m,n],ndir);
g=zeros(data_num,ndir);
% getfig(temps);
for data_numi=1:data_num
    g(data_numi,:)=nansum(nansum(repmat(smp(:,:,data_numi),1,1,ndir).*temps,1),2);
    %     for i=1:ndir
    %         g(data_numi,i)=nansum(nansum(smp(:,:,data_numi).*temps(:,:,i)));
    %     end
end
GC0=cell(1,data_num);
dir2i=0;
dir21i=0;
dir2i1=zeros(data_num,1);
dir21i1=zeros(data_num,1);
for data_numi=1:data_num
    if sum(g(data_numi,:)==0)==size(g,2)
        dir2i1(data_numi)=nan;
        dir21i1(data_numi)=nan;
        continue;
    end
    
    [~,GC0{data_numi}]=max(abs(g(data_numi,:)));
    
    GC=GC0{data_numi};
    if length(GC)>1
        g2=zeros(1,length(GC));
        g2i=0;
        for gc=reshape(GC,1,[])
            g2i=g2i+1;
            tempsi=temps(:,:,gc);
            tempsi(tempsi==-1)=0;
            g2(g2i)=nansum(nansum(smp(:,:,data_numi).*tempsi));
        end
        [~,GC]=max(abs(g2));
        if length(GC)>1
            GC=GC(1);
        end
    end
    dir1=(360/ndir)*(GC-1);
    dir2=dir1+90;if dir2>=360; dir2=dir2-360;end
    dir21=dir1+90+180;if dir21>=360; dir21=dir21-360;end
    if dir2i==0&&dir21i==0
        dir2i=cosd(dir2)+sqrt(-1)*sind(dir2);
        dir21i=cosd(dir21)+sqrt(-1)*sind(dir21);
        dir2i1(data_numi)=dir2;
        dir21i1(data_numi)=dir21;
    else
        if abs(dir2i+cosd(dir2)+sqrt(-1)*sind(dir2))>sqrt(2)
            dir2i=dir2i+cosd(dir2)+sqrt(-1)*sind(dir2);
            dir21i=dir21i+cosd(dir21)+sqrt(-1)*sind(dir21);
            dir2i=dir2i/abs(dir2i);
            dir21i=dir21i/abs(dir21i);
            dir2i1(data_numi)=dir2;
            dir21i1(data_numi)=dir21;
        else
            dir2i=dir2i+cosd(dir21)+sqrt(-1)*sind(dir21);
            dir21i=dir21i+cosd(dir2)+sqrt(-1)*sind(dir2);
            dir2i=dir2i/abs(dir2i);
            dir21i=dir21i/abs(dir21i);
            dir2i1(data_numi)=dir21;
            dir21i1(data_numi)=dir2;
        end
    end
end
dir2i1(isnan(dir2i1))=[];
dir21i1(isnan(dir21i1))=[];
dir2i=mediandir(dir2i1);
dir21i=mediandir(dir21i1);
dir2=dir2i;if dir2<0;dir2=dir2+360;end
dir21=dir21i;if dir21<0;dir21=dir21+360;end

ind2=closeddirind(dir2);
ind21=closeddirind(dir21);
% mod33_temp=zeros(3,3);
% mod33_temp([8,9,6,3,2,1,4,7])=0:45:360-45;
% mod33_temp(2,2)=360;
% ind2=find(abs(dir2-mod33_temp(:))==min(abs(dir2-mod33_temp(:))));
% ind21=find(abs(dir21-mod33_temp(:))==min(abs(dir21-mod33_temp(:))));

mod330=getmod33(smp);
% mod330sum=sum(sum(mod330==0));
% mod330(:,:,mod330sum>=3)=0;
[ind2r,ind2c]=ind2sub([3,3],sort(ind2,'descend'));
[ind21r,ind21c]=ind2sub([3,3],sort(ind21,'ascend'));
if length(ind2r)>1&&length(ind2r)==length(ind21r)
    d221r=zeros(length(ind2r),2);
    for i=1:length(ind2r)
        ind2ri=ind2r(i);
        ind2ci=ind2c(i);
        ind21ri=ind21r(i);
        ind21ci=ind21c(i);
        ind2vi=mod330(ind2ri,ind2ci,:);if sum(ind2vi~=0)==0;ind2vi(:)=99999999;end
        ind21vi=mod330(ind21ri,ind21ci,:);if sum(ind2vi~=0)==0;ind21vi(:)=99999999;end
        ind22vi=mod330(2,2,:);
        d221r(i,:)=[sum((ind2vi-ind22vi).^2),sum((ind21vi-ind22vi).^2)];
    end
    [im,~]=find(d221r==min(d221r(:)));
    im=unique(im);
    if length(im)>1
        im2=find(d221r(im,:)==max(d221r(im,:)));
        im=im(im2(1));
    end
    ind2r=ind2r(im);
    ind2c=ind2c(im);
    ind21r=ind21r(im);
    ind21c=ind21c(im);
else
    ind2r=ind2r(1);
    ind2c=ind2c(1);
    ind21r=ind21r(1);
    ind21c=ind21c(1);
end
mod330(mod330==0)=nan;
ind22v=nanmedian(mod330(2,2,:),3);
ind2v=nanmedian(mod330(ind2r,ind2c,:),3);
ind21v=nanmedian(mod330(ind21r,ind21c,:),3);

if ~isnan(ind22v)
    mod330(:,:,isnan(mod330(2,2,:)))=nan;
else
    dw=0;
    while (isnan(ind22v)||isnan(ind2v)||isnan(ind21v))&&dw<ceil(size(smp,1)/4)
        dw=dw+1;
        mod330=getmod33(smp,dw);
        mod330(mod330==0)=nan;
        ind22v=nanmedian(mod330(2,2,:),3);
        ind2v=nanmedian(mod330(ind2r,ind2c,:),3);
        ind21v=nanmedian(mod330(ind21r,ind21c,:),3);
    end
    while (isnan(ind2v)||isnan(ind21v))&&dw<ceil(size(smp,1)/4)
        dw=dw+1;
        mod330=getmod33(smp,dw);
        mod330(mod330==0)=nan;
        ind2v=nanmedian(mod330(ind2r,ind2c,:),3);
        ind21v=nanmedian(mod330(ind21r,ind21c,:),3);
    end
    if isnan(ind2v);ind2v=99999999;end
    if isnan(ind21v);ind21v=99999999;end
end

if sum((ind2v-ind22v).^2)<=sum((ind21v-ind22v).^2)
    dir11=dir2+90;
    if dir11>=360
        dir11=dir11-360;
    end
    GC=round(dir11/360*ndir+1);
else
    dir11=dir2-90;
    if dir11<0
        dir11=dir11+360;
    end
    GC=round(dir11/360*ndir+1);
end
if GC>ndir
    GC=1;
end
if GC<1
    GC=1;
end
tempsi=temps(:,:,GC);
tempsi(tempsi==-1)=0;
smp(isnan(smp))=0;
smp1=smp.*repmat(tempsi,1,1,data_num);
end
function mod33=getmod33(smp,dw)
if nargin==1
    dw=0;
end
smp(smp==0)=nan;
[m,n,data_num]=size(smp);
mod33=zeros(3,3,data_num);
% olpw2=0;
i1=(1:3)';
ii1=repmat(round((i1-1)/3*m+1),1,3);
ii2=repmat(round(i1/3*m+1),1,3)-1;
j1=(1:3);
jj1=repmat(round((j1-1)/3*n+1),3,1);
jj2=repmat(round(j1/3*n+1),3,1)-1;
% ii1=ii1-olpw2;
% ii2=ii2+olpw2;
% jj1=jj1-olpw2;
% jj2=jj2+olpw2;
ii1=ii1-dw;
ii2=ii2+dw;
jj1=jj1-dw;
jj2=jj2+dw;
ii1(ii1<1)=1;
ii2(ii2>m)=m;
jj1(jj1<1)=1;
jj2(jj2>n)=n;

for k=1:data_num
    mk=nanmedian(nanmedian(smp(:,:,k)));
    for i=1:3
        for j=1:3
            ii=ii1(i,j):ii2(i,j);
            jj=jj1(i,j):jj2(i,j);
            smpii=smp(ii,jj,k);
            %             mkij=wls(smpii(~isnan(smpii)),ones(length(smpii(~isnan(smpii))),2));
            mkij=nanmedian(smpii(:));
            mod33(i,j,k)=mkij/mk;
            if isnan(mod33(i,j,k))
                mod33(i,j,k)=0;
            end
        end
    end
%     if sum(sum(mod33(:,:,k)==0))>=3
%         mod33(:,:,k)=zeros(3);
%     end
end
end
function temps=gettemps(siz,ndir)
temps=zeros(siz(1),siz(2),ndir);
r0=round(siz(1)/2);
c0=round(siz(2)/2);
r1=siz(1)-r0;%the bottom
r2=-r0+1;%the top
c1=siz(2)-c0;%the right
c2=-c0+1;%the left
[x,y]=meshgrid(c2:1:c1,r2:1:r1);
dirs=(0:ndir-1)*(360/ndir);
i=0;
for dir=dirs
    i=i+1;
    [x1,y1,x2,y2]=getintpoint(dir,r1,r2,c1,c2);
    xx1=x1-x;
    yy1=y1-y;
    xx2=x2-x;
    yy2=y2-y;
    xxyy1=[xx1(:),yy1(:),zeros(length(xx1(:)),1)];
    xxyy2=[xx2(:),yy2(:),zeros(length(xx2(:)),1)];
    cro=cross(xxyy1,xxyy2);
    tempsi=zeros(siz);
    tempsi(cro(:,3)>=0)=1;
    tempsi(tempsi==0)=-1;
    temps(:,:,i)=tempsi;
end
end

function [x1,y1,x2,y2]=getintpoint(dir,r1,r2,c1,c2)
%the slope
k=tand(dir);

if k==inf
    k=9999;
end
if k==-inf
    k=-9999;
end

if (k*c1)<=r1&&(k*c1)>=r2 % the interact point is in the right edge
    x1=c1;
    y1=k*c1;
end
if (k*c1)<r2 % the interact point is in the right bottom edge
    x1=r2/k;
    y1=r2;
end
if (k*c1)>r1 % the interact point is in the right top edge
    x1=r1/k;
    y1=r1;
end

if (k*c2)<=r1&&(k*c2)>=r2 % the interact point is in the left edge
    x2=c2;
    y2=k*c2;
end
if (k*c2)<r2 % the interact point is in the left bottom edge
    x2=r2/k;
    y2=r2;
end
if (k*c2)>r1 % the interact point is in the left top edge
    x2=r1/k;
    y2=r1;
end
if dir==90
    x1=0;y1=r1;
    x2=0;y2=r2;
end
if dir==270
    x1=0;y1=r2;
    x2=0;y2=r1;
end
%
dir1=atan2(y1-y2,x1-x2)/pi*180;
if dir1<0
    dir1=dir1+360;
end
if abs(dir1-dir)>1e-6
    x11=x1;y11=y1;
    x1=x2;y1=y2;
    x2=x11;y2=y11;
end
end

function [unk,Llnorm]= wls(L,B,nom)
%wls1
%
%--this code is used to estimate the unknow in L1-norm
%
%Input:
%  L:the observations
%  B:the design matrix
%  nom: 1 L1-norm/2 L2-norm
%Output:
%  unk:the estimated unknowns
%  Llnorm:the sum of residuals
if nargin==2
    nom=2;
end
if nom==1
    n = length(L);
    p = size(B,2);
    f= [zeros(2*(p),1);ones(2*n,1)] ;
    b = L;
    A = [B , -B , eye(n), - eye(n)];
    lb = zeros(2*(n+p),1);
    options = optimoptions('linprog','Display','none','MaxIter', 10000);
    %options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07,'MaxIter', 100);
    [x,fval]= linprog(f,[],[],A,b,lb,[],[],options);
    if length(x)<2*p
        unk=zeros(p,1);
        Llnorm=0;
    else
        unk = x(1:p)- x(p+1:2*(p));
        Llnorm = sum (abs(L-B*unk));
        %         if abs(fval-Llnorm )>1.0e-6
        %             [x,fval]= linprog(f,[], [], A , b, lb, [] , [], options);
        %             unk = x(1:p)-x(p+1:2*(p));
        %             Llnorm = sum(abs(L-B*unk)) ;
        %         end
    end
else
    unk=pinv(B'*B)*B'*L;
    Llnorm = sum(abs(L-B*unk)) ;
end
end
function dir=mediandir(dirs)
maxdir=max(dirs);
mindir=min(dirs);
if maxdir-mindir>180
    for i=1:length(dirs)
        if maxdir-dirs(i)>180
            dirs(i)=dirs(i)+360;
        end
    end
end
dir=median(dirs);
if dir>=360
    dir=dir-360;
end
end
function ind2=closeddirind(dir2)

mod33_temp=zeros(3,3);
mod33_temp([8,9,6,3,2,1,4,7])=0:45:360-45;
mod33_temp(2,2)=500;
if dir2>180
    mod33_temp(8)=360;
end
ind2=find(abs(dir2-mod33_temp(:))==min(abs(dir2-mod33_temp(:))));
end
