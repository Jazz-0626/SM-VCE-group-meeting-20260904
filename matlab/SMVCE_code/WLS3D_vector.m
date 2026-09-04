function Result_wls=WLS3D(pDATA)
%WLS_3D
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to solve 3-D deformations based on WLS method
%
%Inputs:
%  pDATA: the data file including the directory
%       data:       the row*col*data_num matrix, indcates the input SAR displacement
%       inc,azi:    same as data, but the data's incident and azimuth angle
%       losazienu:  1*data_num vector, indicates the corresponding data is los azi or enu data
%       leftorright:1*data_num vector, indicates the corresponding data is right- (1) or left-looking (-1) mode
%       pmask:      the path of mask data, mask is a row*col matrix, 0 means the pixel would not be solved 3-D disp.
%
%Output:
%  Result:a structure variable
%     enu:          the 3-D deformations
%     var:          the variance of observations and 3-D deformations
%     obsP:         the weight of observations determined based on a 5*5 window
%     total_time:   the used time%
%
% Last modified: Sept. 20th, 2021
fprintf('start solving the 3-D deformations through WLS...\n');
start_time=clock;
load(pDATA,'data','inc','azi','losazienu',...
    'leftorright','mask');
ft=0;
if ft==1
    load maskt
    mask=maskt;
end

load('enangle1s.mat','enangle1');
enangle=enangle1;
[row,col,data_num]=size(data);
enanglesin=sin(enangle);
windowsize=15;

fprintf('start calculating the weight of each obs. based on the hypothesis of ergodicity...\n');
P=getWLSP(data,windowsize);
% P=ones(row,col,data_num);
data(isnan(data))=0;
defo_e=zeros(row,col);
defo_n=zeros(row,col);
defo_u=zeros(row,col);
var_e=zeros(row,col);
var_n=zeros(row,col);
var_u=zeros(row,col);
var_obs=zeros(row,col,data_num);

Bgeo=lk_vec(azi,inc,losazienu,leftorright);

mask0=sum(mask,2);
i_start=min(find(mask0~=0));
clear mask0
flag=matlabversion();
if flag==1
    progressBar= CommandLineProgressBar(row);
    progressBar.message = 'Solving 3-D dsip. by WLS: ';
    progressBar.barLength = 42;
else
    progressBar=0;
end
mask=double(mask);
parfor i=1:row
    if flag==1
        progressBar.increment;
    end
    tic
    for j=1:col
        if mask(i,j)==0
            continue;
        end
        L_i=reshape(data(i,j,:),[],1);
        Pi=reshape(P(i,j,:),[],1);
        w2=30;
        ii=max(1,i-w2):min(i+w2,row);
        jj=max(1,j-w2):min(j+w2,col);
        Bgeo_i=reshape(nanmean(nanmean(Bgeo(ii,jj,:))),3,[])';
        
        L_i=[L_i;1e-6];
        tanv=tan(enangle(i,j));
        Bv=[tanv,-1,0];
        tanvs=3;
        if abs(tanv)>tanvs
            Bv=[tanv/abs(tanv)*tanvs,-1,0];
        end
        Bgeo_i=[Bgeo_i;Bv];
        Pi=[Pi;1];
        
        %         if sum(abs(L_i([3,6,9])))==0
        %             Lii=data(ii,jj,[3,6,9]);
        %
        %             Lii(Lii==0)=nan;
        %             L_i([3,6,9])=reshape(nanmean(nanmean(Lii)),3,1);
        %             L_i(isnan(L_i))=0;
        %         end
        
        neq0=find(L_i~=0);
        L_i(end)=0;
        L_i=reshape(L_i(neq0),[],1);
        Bgeo_i=Bgeo_i(neq0,:);
        Pi=Pi(neq0);
        P_i=diag(Pi);
        if rank(Bgeo_i)<3
            continue;
        end
        NN=inv(Bgeo_i'*P_i*Bgeo_i);
        x=NN*Bgeo_i'*P_i*L_i;
        
        if sign(x(2)*enanglesin(i,j))<0
            L_i=reshape(data(i,j,:),[],1);
            Pi=reshape(P(i,j,:),[],1);
            w2=30;
            ii=max(1,i-w2):min(i+w2,row);
            jj=max(1,j-w2):min(j+w2,col);
            Bgeo_i=reshape(nanmean(nanmean(Bgeo(ii,jj,:))),3,[])';
            
            L_i=[L_i;1e-6];
            tanv=-tan(enangle(i,j));
            Bv=[tanv,-1,0];
            tanvs=3;
            if abs(tanv)>tanvs
                Bv=[tanv/abs(tanv)*tanvs,-1,0];
            end
            Bgeo_i=[Bgeo_i;Bv];
            Pi=[Pi;1];
            
            %         if sum(abs(L_i([3,6,9])))==0
            %             Lii=data(ii,jj,[3,6,9]);
            %
            %             Lii(Lii==0)=nan;
            %             L_i([3,6,9])=reshape(nanmean(nanmean(Lii)),3,1);
            %             L_i(isnan(L_i))=0;
            %         end
            
            neq0=find(L_i~=0);
            L_i(end)=0;
            L_i=reshape(L_i(neq0),[],1);
            Bgeo_i=Bgeo_i(neq0,:);
            Pi=Pi(neq0);
            P_i=diag(Pi);
            if rank(Bgeo_i)<3
                continue;
            end
            NN=inv(Bgeo_i'*P_i*Bgeo_i);
            x=NN*Bgeo_i'*P_i*L_i;
            
        end
        
        defo_e(i,j)=x(1);
        defo_n(i,j)=x(2);
        defo_u(i,j)=x(3);
        
        var_e(i,j)=NN(1);
        var_n(i,j)=NN(5);
        var_u(i,j)=NN(9);
    end
    if flag==0
        fprintf(['row:' num2str(i) '/' num2str(row) '; used time for this row:' num2str(toc) 's\n']);
    end
end
if ft==1
    load Result_WLS_20220120120137
    for i=1:3
        di=Result_wls.enu(:,:,i);
        if i==1
            defo_e(mask==0)=di(mask==0);
        elseif i==2
            defo_n(mask==0)=di(mask==0);
        else
            defo_u(mask==0)=di(mask==0);
        end
    end
end

InputData=load(pDATA,'data','inc','azi','losazienu',...
    'leftorright','mask');
total_time=etime(clock,start_time);
Result_wls.enu=cat(3,defo_e,defo_n,defo_u);
Result_wls.var.obs=var_obs;
Result_wls.var.enu=cat(3,var_e,var_n,var_u);
Result_wls.obsP=P;
Result_wls.InputData=InputData;
Result_wls.total_time=total_time;
Result_wls.enu(Result_wls.enu==0)=nan;
Result_wls.var.obs(Result_wls.var.obs==0)=nan;
Result_wls.var.enu(Result_wls.var.enu==0)=nan;
if ft==0
    save(['Result_WLS_',datestr(clock,'yyyymmddHHMMSS')],'Result_wls','-v7.3');
end
fprintf('solving the 3-D deformations through WLS has done...\n');
fprintf('total used time: %f s...\n\n',total_time);
end
%% other functions
function P=getWLSP(data,windowsize)
%getWLSP
%--By LiuJH 2019/03/23
%--This function is used to get the weiht of data pixel-by-pixel on the
%       hypothesis of ergodicity
%
%Inputs:
%  data:        the used data with size of m*n*data_num
%  windowsize:  the window size used for estimating variance of different
%               observations
%Output:
%  P   : the weiht of data used for WLS
%
% Last modified: Sept. 20th, 2019

data(data==0)=nan;
[row,col,datanum]=size(data);
y=zeros(row,col,datanum);
windowsize2=round((windowsize-1)/2);
[xCol,yRow]=meshgrid(1:col,1:row);

[x1,y1]=meshgrid(-windowsize2:windowsize2,-windowsize2:windowsize2);
temp1=repmat(x1(:)',row*col,1);
temp2=repmat(y1(:)',row*col,1);
xCol1=repmat(xCol(:),1,((windowsize2*2)+1)^2)+temp1;
yRow1=repmat(yRow(:),1,((windowsize2*2)+1)^2)+temp2;
flag=(xCol1>0)+(xCol1<=col)+(yRow1>0)+(yRow1<=row);

ind=zeros(row*col,((windowsize2*2)+1)^2);
ind(flag==4)=sub2ind([row,col],yRow1(flag==4),xCol1(flag==4));
clear temp1 temp2 yRow1 xCol1 xCol yRow
for i=1:datanum
    datatemp=squeeze(data(:,:,i));
    datatempi=zeros(row*col,((windowsize2*2)+1)^2);
    datatempi(ind>0)=datatemp(ind(ind>0));
    datatempi(datatempi==0)=nan;
    y(:,:,i)=reshape(nanstd(datatempi'),row,col);
end
P=1./y;
P(y==0)=0;
P(isnan(y))=nan;
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
