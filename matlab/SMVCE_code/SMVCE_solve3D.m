function Result_smvce=SMVCE_solve3D(pDATA)
%SMVCE_solve3D
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%--This function is used to solve 3-D deformations based on SMVCE method
%
%Inputs:
%  pDATA: the data file including the directory
%       data:       the row*col*data_num matrix, indcates the input SAR displacement
%       inc,azi:    same as data, but the data's incident and azimuth angle
%       losazienu:  1*data_num vector, indicates the corresponding data is los azi or enu data
%       leftorright:1*data_num vector, indicates the corresponding data is right- (1) or left-looking (-1) mode
%       dem:        dem data, dem is a row*col matrix
%       mask:       mask data, mask is a row*col matrix, 0 means the pixel would not be solved 3-D disp.
%       coor：      the coordinate of the data
%       windowsize: the window size for strain model
%       fault:     the trace of known fault rupture
%       fsmpara:    the number of the strain model dimension, 2 or 3, default 2
%       flag_BOI:   1*data_num vector, indicates the corresponding data is BOI data
%       flag_smad:  if use the SMAD method to adptively remove the effect in the fault rupture area
%       flag_adpws: if adaptive enlarge the windowsize to include more pixels for the strain model
%       flag_interWeight: if determine the internal weight between pixels of the same kind of SAR disp
%
%Output:
%  Result:a structure variable
%     enu:          the 3-D deformations
%     var:          the variance of observations and 3-D deformations
%     para_sm:      the estimated strain parameters
%     sita:the      convergence result for each pixel
%     SHPcount:     the number of select homogeneous points (SHP) for each pixel
%     total_time:   the used time%
%
% Last modified: Sept. 20th, 2021

DATA=load(pDATA,'data','inc','azi','losazienu',...
    'leftorright','dem','mask','coor',...
    'flag_smad','flag_adpws','flag_interWeight',...
    'windowsize','fault','fsmpara','flag_if_2D');
dem=DATA.dem;
mask=DATA.mask;
flag_smad=DATA.flag_smad;
flag_adpws=DATA.flag_adpws;
flag_interWeight=DATA.flag_interWeight;
fsmpara=DATA.fsmpara;
inc=DATA.inc;
azi=DATA.azi;
losazienu=DATA.losazienu;
leftorright=DATA.leftorright;
data=DATA.data;
flag_BOI=zeros(size(data,3),1);
windowsize0=DATA.windowsize;
coor=DATA.coor;
fault=DATA.fault;
flag_if_2D=DATA.flag_if_2D;
clear DATA

[row,col,data_num]=size(data);

maskfault=zeros(size(mask));
if iscell(fault)
    [cs,rs]=meshgrid(1:col,1:row);
    for i=1:length(fault)
        fi=lonlat2sub(fault{i},coor);
        plyb=polybuffer(fi(:,2:-1:1),'lines',windowsize0*2.5);
        polyi=[plyb.Vertices;plyb.Vertices(1,:)];
        [in,on]=inpolygon(cs,rs,polyi(:,1),polyi(:,2));
        maskfault=maskfault+in+on;
    end
end
fault0=fault;


data=double(data);
data(isnan(azi))=nan;
mask=double(mask);
dem=double(dem);

data(isnan(data))=0;
% data1=data;
% data1(data1==0)=nan;
% ax1=getfig(data1,datalim(data1));
% for j=1:size(data1,3)
%     axes(ax1(j));
%     if length(fault(:))~=1
%     for i=1:length(fault)
%         subi=lonlat2sub(fault{i},coor);
%         hold on
%         plot(subi(:,2),subi(:,1));
%     end
%     end
% end


unknum=12;
if fsmpara==2&&flag_if_2D==0
    unknum=9;
    npara_sm=6;
end
if fsmpara==2&&flag_if_2D==1
    unknum=6;
    npara_sm=4;
end
if fsmpara==3&&flag_if_2D==0
    unknum=12;
    npara_sm=9;
end
if fsmpara==3&&flag_if_2D==1
    unknum=8;
    npara_sm=6;
end


defo_e=zeros(row,col);
defo_n=zeros(row,col);
defo_u=zeros(row,col);
para_sm=zeros(row,col,npara_sm);
var_e=zeros(row,col);
var_n=zeros(row,col);
var_u=zeros(row,col);
var_obs=zeros(row,col,data_num);
[xx,yy]=meshgrid(1:col,1:row);
xy_m_s=1000;
x_m=abs(coor.post_lon)*108*xy_m_s;
y_m=abs(coor.post_lat)*108*xy_m_s;
sitaxy=zeros(row,col,data_num);
SHPcount=ones(row,col,data_num);

if exist('sm_vce.log','file')
    delete('sm_vce.log');
end
Bgeo=lk_vec(azi,inc,losazienu,leftorright);
fprintf('start solving the 3-D deformations through SM-VCE...\n');

start_time=clock;
flag=matlabversion();
% flag=0;
if flag==1
progressBar= CommandLineProgressBar(row);
progressBar.message = 'Solving 3-D dsip. by SM-VCE: ';
progressBar.barLength = 42;
else 
progressBar=0;
end

%the data that with larger windowsize
ord_window=3;
ind_maxwindow=[];
ind_maxwindow0=setdiff(1:size(data,3),ind_maxwindow);
data0=data(:,:,ind_maxwindow0);
data1=data(:,:,ind_maxwindow);
Ts=[];
parfor i=1:row
    if flag==1
    progressBar.increment;
    end
    tic
    for j=1:col
        if mask(i,j)==0
            continue;
        end
        windowsize=windowsize0;

        maxwin=5*windowsize;
        maxdis=(maxwin/2)*sqrt(x_m^2+y_m^2)/1000;
        diswa=1/(1-maxdis);
        diswb=maxdis/(maxdis-1);
        
        %通过断层线fault或者自适应flag_smad去除断层另一侧的点的算法，或者自适应窗口大小变化flag_adpws得到相应的观测数据
        
        if maskfault(i,j)>0
            fault=fault0;
        else
            fault=0;
        end
        t0i=clock;
        if ~isempty(ind_maxwindow)
            [L_i00,iis0,jjs0,SHPcounti0]=getL(data0,i,j,windowsize,coor,fault,flag_smad,flag_adpws);
            [L_i01,iis1,jjs1,SHPcounti1]=getL(data1,i,j,windowsize*ord_window,coor,fault,flag_smad,flag_adpws);
            
            L_i0=cell(1,data_num);
            iis=zeros(data_num,2);
            jjs=zeros(data_num,2);
            SHPcounti=zeros(1,1,data_num);
            L_i0(ind_maxwindow0)=L_i00;
            L_i0(ind_maxwindow)=L_i01;
            
            iis(ind_maxwindow0,:)=iis0;
            iis(ind_maxwindow,:)=iis1;
            jjs(ind_maxwindow0,:)=jjs0;
            jjs(ind_maxwindow,:)=jjs1;
            SHPcounti(ind_maxwindow0)=SHPcounti0;
            SHPcounti(ind_maxwindow)=SHPcounti1;
        else
            [L_i0,iis,jjs,SHPcounti]=getL(data,i,j,windowsize,coor,fault,flag_smad,flag_adpws);
        end
        t0i1=etime(clock,t0i);
        t0i=clock;
        
        SHPcount(i,j,:)=SHPcounti;
        %得到用于方差分分量估计的观测值，设计矩阵和权重矩阵
        
        [L,B,k,P]=getLBkP(L_i0,iis,jjs,Bgeo,xx,yy,x_m,y_m,fsmpara,dem,i,j,flag_interWeight,flag_BOI,diswa,diswb,flag_if_2D);
        t0i2=etime(clock,t0i);
        t0i=clock;

        w2=5;
        ii=max(1,i-w2):min(i+w2,row);
        jj=max(1,j-w2):min(j+w2,col);
        Bgeo_i=reshape(nanmean(nanmean(Bgeo(ii,jj,:))),3,[])';
        neq0=find(k>=0.1*windowsize^2);
        Bgeo_i=Bgeo_i(neq0,:);
        if flag_if_2D==0
            rankB=3;
        else
            rankB=2;
            Bgeo_i(:,2)=[];
        end
        if rank(Bgeo_i)<rankB||cond(Bgeo_i)>1e10
            continue;
        end
%         if sum(k>=0.1*windowsize^2)<3 %||rank(B0'*B0)<unknum
%             continue;
%         end
        
        %方差分量估计
        %VCE DETERMINE THE WEIGHT
        [sita,P_vce,unkn]=SMVCE_vce(L,B,P);
        
        t0i3=etime(clock,t0i);
        t0i=clock;
%         a=cell2mat(reshape(L_i0,1,1,[]));
%         a(a==0)=nan;
%         getfig(a);
        
        sitaxy(i,j,:)=sita;
        %将定权失败的点输出到日志文件
        
        if sum(sita)==0
            fid=fopen('sm_vce.log','a');
            fprintf(fid,'row:%d,col:%d are not calculated for the  divergence of VCE\n',i,j);
            fclose(fid);
        end
        
        %各类观测值的权重
        %the variance of each observation
        p=zeros(data_num,1);
        var_obs_i=zeros(data_num,1);
        for kk=1:data_num
            pii=diag(P_vce{kk});
            pii(pii==0)=nan;
            p(kk)=nanmean(pii);
            if isnan(p(kk))
                p(kk)=0;
                continue;
            end
            var_obs_i(kk)=sita(kk)/(p(kk)+eps);
        end        
        var_obs(i,j,:)=sqrt(var_obs_i(:));
        
        %得到的待求参数
        if flag_if_2D==0
        defo_e(i,j)=unkn(1);
        defo_n(i,j)=unkn(2);
        defo_u(i,j)=unkn(3);
        para_sm(i,j,:)=unkn(4:end);
        
        %三维形变结果的方差
        %the variance of 3-D deformations
        Bgeo_var=reshape(Bgeo(i,j,:),3,data_num)';
        Bgeo_var(isnan(Bgeo_var))=0;
        sita(sita==0)=nan;
        d0=nanmean(sita);
        var_enu=d0*pinv(Bgeo_var'*diag(p)*Bgeo_var);
        var_e(i,j)=sqrt(var_enu(1));
        var_n(i,j)=sqrt(var_enu(5));
        var_u(i,j)=sqrt(var_enu(9));
        else
        defo_e(i,j)=unkn(1);
        defo_u(i,j)=unkn(2);
        para_sm(i,j,:)=unkn(3:end);
        
        %三维形变结果的方差
        %the variance of 3-D deformations
        Bgeo_var=reshape(Bgeo(i,j,:),3,data_num)';
        Bgeo_var(isnan(Bgeo_var))=0;
        sita(sita==0)=nan;
        d0=nanmean(sita);
        Bgeo_var(:,2)=[];
        var_enu=d0*pinv(Bgeo_var'*diag(p)*Bgeo_var);
        var_e(i,j)=sqrt(var_enu(1));
        var_u(i,j)=sqrt(var_enu(4));

        end
        t0i4=etime(clock,t0i);
%         Ts=[Ts;t0i1,t0i2,t0i3,t0i4];
    end
    if flag==0
    fprintf(['row:' num2str(i) '/' num2str(row) '; used time for this row:' num2str(toc) 's\n']);
    end
end
total_time=etime(clock,start_time);

enu=cat(3,defo_e,defo_n,defo_u);
para_sm1=para_sm;
var_enu=cat(3,var_e,var_n,var_u);
var_obs1=var_obs;
sitaxy1=sitaxy;
SHPcount1=SHPcount;

enu(enu==0)=nan;
para_sm1(para_sm1==0)=nan;
var_enu(var_enu==0)=nan;
var_obs1(var_obs1==0)=nan;
sitaxy1(sitaxy1==0)=nan;
SHPcount1(SHPcount1==0)=nan;

InputData=load(pDATA);
Result_smvce.enu=enu;
Result_smvce.var.obs=var_obs1;
Result_smvce.var.enu=var_enu;
Result_smvce.para_sm=para_sm1;
Result_smvce.sita=sitaxy1;
Result_smvce.coor=coor;
Result_smvce.InputData=InputData;
Result_smvce.SHPcount=SHPcount1;
Result_smvce.total_time=total_time;
% save(['Result_SMVCE_',datestr(clock,'yyyymmddHHMMSS')],'Result_smvce','-v7.3');
total_time=etime(clock,start_time);
fprintf('solving the 3-D deformations through SM-VCE has done...\n');
fprintf('total used time: %f s...\n',total_time);
end


%% other functions
function [L,B,k,P]=getLBkP(L_i0,iis,jjs,Bgeo,xx,yy,x_m,y_m,fsmpara,dem,i,j,flag_interWeight,flag_BOI,diswa,diswb,flag_if_2D)
data_num=length(L_i0);
L=cell(data_num,1);
B=cell(data_num,1);
k=zeros(data_num,1);
P=cell(data_num,1);
for kk=1:data_num
    iii=iis(kk,1):iis(kk,2);
    jjj=jjs(kk,1):jjs(kk,2);
    Bgeo_i=Bgeo(iii,jjj,kk*3-2:kk*3);
    dx_i=x_m*(xx(iii,jjj)-j);
    dy_i=y_m*(yy(iii,jjj)-i);
    if fsmpara==3
        dz_i=(dem(iii,jjj)-dem(i,j));
        if flag_if_2D==0
        B_i=get_design_mat(Bgeo_i,dx_i,dy_i,dz_i);
        else

        B_i=get_design_mat_2D(Bgeo_i,dx_i,dy_i,dz_i);
        end
    else
        if flag_if_2D==0
        B_i=get_design_mat(Bgeo_i,dx_i,dy_i);
        else
        B_i=get_design_mat_2D(Bgeo_i,dx_i,dy_i);
        end
    end
    tem=L_i0{kk};
    tem_B=B_i;        
    neq0=find(tem(:).*(~isnan(sum(B_i,2)))~=0);
    L{kk}=tem(neq0);
    B{kk}=tem_B(neq0,:);
    k(kk)=length(neq0);
    if flag_interWeight==1
        Li=L{kk};
        Bi=[dx_i(neq0)/1000,dy_i(neq0)/1000,ones(length(neq0),1)];
        if k(kk)<1.5*length(iii)
            wi=ones(1,k(kk));
        else
            [~,wi]=iwlsnew(Bi,Li,ones(length(neq0),1));
        end
    else
        wi=ones(1,k(kk));
    end
    if flag_BOI(kk)==1
        disf=max(1,sqrt(dx_i.^2+dy_i.^2)/1000);
        wi=wi.*reshape(diswa*disf(neq0)+diswb,size(wi));
    end
    P{kk}=sparse(1:k(kk),1:k(kk),wi,k(kk),k(kk));
end
end
function [x1,wi]=iwlsnew(B,L,wi,flagwls1,ita,itanum,flsmr)
% tt=zeros(1,2);
% tici=0;
% tic
t1=1.5;
t2=2.5;
if ~exist('flagwls1','var')
    flagwls1=0;
end
if ~exist('itanum','var')
    itanum=10;
end
if ~exist('flsmr','var')
    flsmr=0;
end
if ~exist('ita','var')
    ita=1e-1;
end
if flagwls1
    [~,wiscale]=wls1s(B,L);
    wi=wi.*wiscale;
end
c=4.685;
nL=length(L);
nx=size(B,2);

W=sparse(1:nL,1:nL,wi,nL,nL);
BW=B'*W;
x=(BW*B)\BW*L;
v=B*x-L;
x1=x;
x(:)=10000;
k=0;
c=median(abs(v));
while sum(abs(x1-x)>ita)>0
    k=k+1;
    x=x1;
    if k>itanum
        break;
    end
    v=B*x-L;
    g=ones(size(v));
    g(abs(v)>t2*c)=0;
    ind=find((abs(v)>t1*c).*(abs(v)<=t2*c)==1);
    g(ind)=1./(abs(v(ind))+1-t1*c);
    wi=wi.*g;
    W=sparse(1:nL,1:nL,wi,nL,nL);
    BW=B'*W;
    BWB=(BW*B);
    BWB(isnan(BWB))=0;
    BWB(abs(BWB)>10^16)=0;
    if flsmr==0
        if rank(full(BWB))<size(BWB,2)
            x1=lsmr(BWB+eye(size(BWB,2))*10^-6,BW*L);
        else
            x1=(BWB\BW*L);
        end
    else
        x1=lsmr(BWB,BW*L);
    end
end
end
function [x1,wi]=iwls(B,L,wi,flagwls1,ita,itanum,flsmr)
nL1=length(L);
ind1=randperm(length(L),nL1);
if ~exist('flagwls1','var')
    flagwls1=0;
end
if ~exist('itanum','var')
    itanum=3;
end
if ~exist('flsmr','var')
    flsmr=0;
end
if ~exist('ita','var')
    ita=1e-1;
end
if flagwls1
    [~,wiscale]=wls1s(B,L);
    wi=wi.*wiscale;
end
c=4.685;
nL=length(L);
nx=size(B,2);

W=sparse(1:nL1,1:nL1,wi(ind1),nL1,nL1);
BW=B(ind1,:)'*W;
x=(BW*B(ind1,:))\BW*L(ind1,:);
% iBB=speye(nx)/(B'*B);
% h=sum((B*iBB).*B,2);
x1=x;
x(:)=10000;
k=0;
while sum(abs(x1-x)>ita)>0
    k=k+1;
    x=x1;
    if k>itanum
        break;
    end
    v=B*x-L;
    s=median(abs(v(wi~=0)));
    %     iBWB=speye(nx)/(B'*B);
    %         h=sum((B*iBWB).*B,2);
%     u=v/s./sqrt(1-h);
    u=v/s;
    g=(1-(u/c).^2).^2;
    g(abs(u)>c)=0;
    wi=wi.*g;
    W=sparse(1:nL1,1:nL1,wi(ind1,:),nL1,nL1);
    BW=B(ind1,:)'*W;
    BWB=(BW*B(ind1,:));
    BWB(isnan(BWB))=0;
    BWB(abs(BWB)>10^16)=0;
    if flsmr==0
        if rank(full(BWB))<size(BWB,2)
            x1=lsmr(BWB+eye(size(BWB,2))*10^-6,BW*L);
        else
            x1=(BWB\BW*L(ind1,:));
        end
    else
        x1=lsmr(BWB,BW*L(ind1,:));
    end
end
end
function [x1,flag0]=wls1s(B,L)
s=3;
x1=B\L;
res=B*x1-L;
medres=median(res);
dres=abs(res-medres);
dres1=sort(dres);
ind=find(dres<dres1(round(0.65*length(L))));
B2=B(ind,:);
L2=L(ind,:);
x2=B2\L2;
res2=B2*x2-L2;
sig0=sqrt(mean(res2.^2));
x1=x2;
flag0=ones(length(L),1);
flag0(abs(res/sig0)>s)=0;
end

function [L_i0,iis,jjs,SHPcounti]=getL(data,i,j,windowsize,coor,fault,flag_smad,flag_adpws)
[row,col,data_num]=size(data);
L_i0=cell(1,data_num);
maxwin=5*windowsize;
minnum=0.2*windowsize*windowsize;
corner_lon=coor.corner_lon;
corner_lat=coor.corner_lat;
post_lon=coor.post_lon;
post_lat=coor.post_lat;

ii=max([1,i-round((windowsize-1)/2)]):min([row,i+round((windowsize-1)/2)]);
jj=max([1,j-round((windowsize-1)/2)]):min([col,j+round((windowsize-1)/2)]);
corner_loni=corner_lon+(min(jj)-1)*post_lon;
corner_lati=corner_lat+(min(ii)-1)*post_lat;
r0=round(median(ii)-ii(1)+1);
c0=round(median(jj)-jj(1)+1);
mask_all_ij=data(i,j,:)~=0;

%different obs may use different windowsize
iis=repmat([min(ii),max(ii)],data_num,1);
jjs=repmat([min(jj),max(jj)],data_num,1);
%基于断层线数据进行同质点选取
%Select the SHPs
if iscell(fault)
    iffault=1;
    mask_fault=getHomoPoints(fault,[length(ii),length(jj)],...
        corner_loni,corner_lati,post_lon,post_lat,r0,c0);
else
    iffault=0;
    mask_fault=ones(length(ii),length(jj));
end
mask_fault=repmat(mask_fault,1,1,data_num);
L_i=data(ii,jj,:).*mask_fault;

%利用SMAD方法自适应选取

if flag_smad==1
    [~,L_i]=SMVCE_SMAD(L_i,r0,c0,1,iffault*(sum(mask_fault(:)==1)~=length(mask_fault(:))));
end

mi_sum=sum(sum(L_i~=0));
mi_sum00=mi_sum;
for L_i0_i=1:data_num
    L_i0{L_i0_i}=L_i(:,:,L_i0_i);
end
windowsizeij=windowsize;
mask_all_ij=mi_sum>0;
%在设定最小观测值个数的前提下，是否需要不断扩大窗口来增加点的个数
%if the correponding obs is not included but must included by
%predefined, then increase the windowsize
if flag_adpws==1
    while sum(mi_sum(mask_all_ij==1)>minnum)~=sum(mask_all_ij)
        mi_sum0=mi_sum;
        windowsizeij=windowsizeij*1.2;
        if windowsizeij>maxwin
            break;
        end
        ii=max([1,i-round((windowsizeij-1)/2)]):min([row,i+round((windowsizeij-1)/2)]);
        jj=max([1,j-round((windowsizeij-1)/2)]):min([col,j+round((windowsizeij-1)/2)]);
        corner_loni=corner_lon+(min(jj)-1)*post_lon;
        corner_lati=corner_lat+(min(ii)-1)*post_lat;
        r0=round(median(ii)-ii(1)+1);
        c0=round(median(jj)-jj(1)+1);
        if iscell(fault)
            mask_fault=getHomoPoints(fault,[length(ii),length(jj)],...
                corner_loni,corner_lati,post_lon,post_lat,r0,c0);
        else
            mask_fault=ones(length(ii),length(jj));
        end
        mask_fault=repmat(mask_fault,1,1,data_num);
        L_i=data(ii,jj,:).*mask_fault;
        if flag_smad==1
            [~,L_i]=SMVCE_SMAD(L_i,r0,c0,1,iffault*(sum(mask_fault(:)==1)~=length(mask_fault(:))));
        end
        mi_sum=sum(sum(L_i~=0));
        idxtemp1=find(mask_all_ij==1);
        idxtemp2=idxtemp1((mi_sum0(idxtemp1)<=minnum).*(mi_sum(idxtemp1)>minnum)==1);
        if ~isempty(idxtemp2)
            iis(idxtemp2,:)=repmat([min(ii),max(ii)],length(idxtemp2),1);
            jjs(idxtemp2,:)=repmat([min(jj),max(jj)],length(idxtemp2),1);
            for L_i0_i=idxtemp2(:)'
                L_i0{L_i0_i}=L_i(:,:,L_i0_i);
            end
            mi_sum00(idxtemp2)=mi_sum(idxtemp2);
        end
    end
end
SHPcounti=mi_sum00;
end
function data=getHomoPoints(faulttraces,siz,corner_lon,corner_lat,post_lon,post_lat,r0,c0)
%getHomoPoints
%--By LiuJH 2019/03/23
%--This function is used to get HomoPoints by fault seperation
%
%Inputs:
%  faulttraces: the fault locations (lon,lat)
%  siz:     the size of the predefined window
%  corner_lon,corner_lat,post_lon,post_lat:the lon,lat of the left/top
%           corner, and the pixel increment
%  r0,c0:   the intersted point
%Output:
%  data:    a 0/1 matrix with size of siz, the pixel with 1 will be further
%           used, and 0 will be obandent
%
% Last modified: June 2019

m=siz(1);
n=siz(2);
data=ones(m,n);

if nargin==6
r0=round((m+1)/2);
c0=round((n+1)/2);
end
% data(r0,c0)=2;
[x,y]=meshgrid(0:n-1,0:m-1);
lons=corner_lon+x*post_lon;
lats=corner_lat+y*post_lat;
lonlim=[corner_lon,corner_lon+(n-1)*post_lon];
latlim=[corner_lat+(n-1)*post_lat,corner_lat];

for i=1:length(faulttraces)
    faulti=faulttraces{i};
    for j=2:size(faulttraces{i},1)
        if ((faulti(j,1)>lonlim(2))&&(faulti(j-1,1)>lonlim(2)))||((faulti(j,1)<lonlim(1))&&(faulti(j-1,1)<lonlim(1)))...
                ||((faulti(j,2)>latlim(2))&&(faulti(j-1,2)>latlim(2)))||((faulti(j,2)<latlim(1))&&(faulti(j-1,2)<latlim(1)))
            continue
        else
            
            x1=min(faulti(j-1,1),faulti(j,1));
            x2=max(faulti(j-1,1),faulti(j,1));
            y1=min(faulti(j-1,2),faulti(j,2));
            y2=max(faulti(j-1,2),faulti(j,2));
            
            x11=min(lons,lons(r0,c0));
            x21=max(lons,lons(r0,c0));
            y11=min(lats,lats(r0,c0));
            y21=max(lats,lats(r0,c0));

            
            k2=(faulti(j,2)-faulti(j-1,2))/(faulti(j,1)-faulti(j-1,1));
            k=(lats-lats(r0,c0))./(lons-lons(r0,c0));
%             k(abs(k)==inf)=99999999;
            b2=faulti(j,2)-k2*faulti(j,1);
            b=lats-k.*lons;
            %get the coordinate of the inter-point of two lines
            intlon=-(b-b2)./(k-k2);
            intlon(r0,:)=(lats(r0,:)-b2)/k2;
            
            intlon(:,c0)=lons(r0,c0);
            
            intlat=-(-b2*k+b*k2)./(k-k2);
            intlat(:,c0)=k2*lons(r0,c0)+b2;
            intlat(r0,:)=lats(r0,c0);
            
            temp=(intlon>x1).*(intlon<x2).*(intlat>y1).*(intlat<y2)... %in polygon
                .*(intlon>=x11).*(intlon<=x21).*(intlat>=y11).*(intlat<=y21); %the inter-point is between the central and other pixels
            data(temp==1)=0;
        end
    end
end
data(r0,c0)=1;
end
function fault=readfault(pfault)
%read exist faults data

fid=fopen(pfault,'r');
faultidex=0;

while ~feof(fid)
    line=fgetl(fid);
    if isempty(line)||strcmp(line,' ')
        continue;
    end
    if contains(line,'>')
        if faultidex==0
            faultidex=faultidex+1;
            coor=[];
        else
            fault{faultidex}=coor;
            faultidex=faultidex+1;
            coor=[];
        end
        continue;
    end
    lines=strsplit(line,',');
    coor=[coor;str2double(lines{1}),str2double(lines{2})];
end
if ~isempty(coor)
    fault{faultidex}=coor;
end
fclose(fid);
end

function lim=datalim(data,ord)
if nargin==1
    ord=3;
end
data(isnan(data))=0;
data(abs(data)==Inf)=0;
a=data(data~=0);
m=nanmean(a(:));
sd=nanstd(a(:));
cmin=m-ord*sd;
cmax=m+ord*sd;
lim=[cmin,cmax];
end

function flag=matlabversion()
str=version();
y=2018;
flag=1;
for i=2000:y-1
   stri=num2str(i);
   if ~isempty(strfind(str,stri))
       flag=0;
       break;
   end
end
end
