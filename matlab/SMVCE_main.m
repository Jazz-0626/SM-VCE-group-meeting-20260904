%There are two SMVCE_main.m SMVCE_main_1.m files, both can calculate 3-D
%displacements of the 8 January 2022 Mw6.7 Menyuan earthquake as an
%example. However, SMVCE_main.m uses less SAR observations with high computational efficiency
%and the obtained 3-D displacements has relatively low accuracy compared with the
%SMVCE_main_1.m, which uses more SAR observation with low computational efficiency but high accuracy.
%If you have a high-equipped computer, you can first try SMVCE_main_1.m,
%else SMVCE_main.m is prefered.
%Users can try these two cases to find more information about this code.
%and the modification of these code is also very welcome.
%Contact me: Dr.Jihong Liu, liujihong@csu.edu.cn
%%
clc;close all;clear
addpath([pwd,filesep,'SMVCE_code']);
load cmap_BlueGreenBrown
set(0,'defaultfigureColormap',cmap_BlueGreenBrown);
set(0,'defaultAxesFontName', 'FreeSans');
set(0,'defaultTextFontName', 'FreeSans');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SM-VCE %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Step 1: view SMVCE_DATA/data_information
%Step 2: if there is DEM data, please put the dem.tif in the folder SMVCE_code, 
%           generally dem.tif is not madatary
%        if there is mask data, please put the mask.tif in the folder SMVCE_code, 
%           the mask file with size of row*col, 1 and 0 indicate the pixels are or
%           are not to calculate the 3-D disp. This file can be preset to exclude
%           those pixels with serious decorrelations. Also, if it is no
%           need to mask area, this file mask.tif is not madatary
%        if the surface is ruptured, and the trace's lon/lat are known, a file 
%           fault.xy should be put in the folder SMVCE_code, two columns. 
%           Also, if there is no known faults, this file fault.xy is not madatary
%Step 3: Check the parameters in the following and run this matlab code.
%步骤1：查看文件SMVCE_DATA/data_information内的介绍并准备相应数据
%步骤2：如果有DEM数据，请将dem.tif放到SMVCE_DATA文件夹下，实际情况中，DEM数据往往不需要参与三维形变解算
%      如果存在掩膜文件mask.tif，大小为row*col, 1代表求解相应像素的三维形变，0代表不求解
%           这个文件是为了不解算低相干区域。需要将这个mask.tif文件也放到SMVCE_DATA文件夹下
%           如果不需要进行掩膜的话，则可以忽略这个文件,或者当这个文件不存在时，程序会解算所有观测数据的重叠区域
%      如果已知地表破裂信息，即断裂线M*2的经纬度信息，这主要是针对地震等地表形变导致的地表跳变
%           可以通过手动勾画等方式获取破裂线的经纬度信息，则需要将fault.xy的文本数据放到SMVCE_DATA文件夹下
%           如果不知道地表破裂信息，则可以忽略这个文件
%步骤3：检查下面的输入参数，运行此代码。
%% calc. the 3-D disp. based on the strain model and variance component estimation (SM-VCE)
%基于地应变模型和方差分量估计获取三维地表形变 

%Put all required data file in the SMVCE_DATA folder
%将所有需要的数据放入文件夹SMVCE_DATA中，包括形变数据的地理tif格式
%以及包含形变观测数据基本信息的data_information文件,以及
%dem.tif(DEM数据，此文件非必需)和mask.tif(掩膜数据，1代表解算，0代表不解算相应点形变，此文件非必需)
%fault.xy代表已知断层线的经纬度信息，两列数据，此文件非必需

pSMVCE_DATA=[pwd,filesep,'SMVCE_DATA'];

%This SMVCE_readdispdata function is to read the Disp measurements file
%读取形变观测数据
%data:          the row*col*N disp. measurements; 形变观测数据
%inc,azi:       N*1 vector, the incident and azimuth angle,入射角和方位角，unit:degree,单位为度
%losazienu:     N*1 vector, the flag indicate the observing geometry of the
%                   cooresponding measurements,1-5 represent LOS, AZI, E-W,
%                   N-S, and Vertical displacements; 代表观测数据的成像几何，1-5
%                   分别代表LOS向，方位向，东西向，南北向和垂直向形变观测数据
%leftorright:   N*1 vector, indicate the measurements are (1) right- or
%                   (-1) left-looking mode data, 代表观测值是右视（1）还是左视（-1）
%coor：          a struct, includes the coordinate information of the data;
%                   包含了数据的坐标信息，例如左上角坐标，及像素坐标增量
%datainfo：      The same information in the file data_information;
%                   和data_information文件中包含的信息是一致的
%dem:           The DEM file, which is not madatary and can be a random matrix
%                   研究区域的DEM数据，非必须，可设为大小为row*col的随机矩阵                    
%mask：         掩膜文件，大小为row*col, 1代表求解相应像素的三维形变，0代表不求解，这个mask文件可以提前
%                   将低相干区域不进行解算,如果不进行掩膜，则将此文件设为row*col的全1矩阵
%                   the mask file with size of row*col, 1 and 0 indicate the pixels are or
%                   are not to calculate the 3-D disp. This file can be preset to exclude
%                   those pixels with serious decorrelations. Also, the
%                   elements can all be 1.
%fault:         the lon,lat of the ruptured fault trace, two columns, which
%                   is not madatary and can be 0；地表破裂线文件，M*2的经纬度信息
%                   主要是针对地震等地表形变导致的地表跳变，通过手动勾画等方式获取破裂线的经纬度信息
%                   此矩阵也可以不设定，直接设成0即可

[data,inc,azi,losazienu,leftorright,coor,dem,mask,fault,datainfo]=SMVCE_readdispdata();

%进行地应变建模的窗口大小
%the window size for modeling the strain model
windowsize=41;

%当只有升降轨LOS向形变时，只能解算东西和垂直向二维形变
%when only ascending and descending LOS observations are available, we can
%only calculate 2D EW/UD displacement
%0 means that we have azimuth observations and we can calculate 3D displacement
%1 means that we don't have azimuth observations and we can calculate 2D displacement
flag_if_2D=0;

%地应变模型的维度，原始情况是三维形变对三维空间求偏导，因此fsmpara=3
%实际应用中三维形变对东西向和南北向空间求偏导即可，所有fsmpara=2
%the dimension of the strain model, originally is 3, but we recomemand 2 in
%real cases
fsmpara=2;

%窗口内同一类观测值之间是否要相对定权 SMVCE_main_2
%if determine the raletive weight between measurements in the window
flag_interWeight=0;

%是否需要满足一定窗口范围内的观测值个数大于某一阈值，而自适应扩大窗口
%if it is neccessary to enlarge the windowsize for more measurements
flag_adpws=0;

%是否需要SMAD算法(环境遥感论文)来自适应剔除断层另一侧的点
%if it is neccessary to exclude the inhomogeneous points in the other side of the fault
flag_smad=0;


indnotuse=[3,6];
dl=[-0.5,0.5];
%show data
datanms=datainfo(:,1);
data(data==0)=nan;
[axs,hcs]=getfig(data,dl,1,1,0,0,datanms);
set(gcf,'Position',[245 56 1484 1232],'name','DATA');
reshapefigure(coor);
for i=1:length(hcs)
    clbtitle(hcs(i),'[m]',14);
    axes(axs(i));hold on
coortick(coor,[0.2,0.2]);
    if ismember(i,indnotuse)
        set(gca,'ycolor','r','xcolor','r');
    end
    if length(fault(:))~=1
        for j=1:length(fault)
            subi=lonlat2sub(fault{j},coor);
            hold on
            plot(subi(:,2),subi(:,1),'m','LineWidth',2);
        end
    end
end
nm='DATA.png';
print(gcf,nm,'-dpng','-r300');
rmfigbg(nm);


data(:,:,indnotuse)=[];
inc(:,:,indnotuse)=[];
azi(:,:,indnotuse)=[];
losazienu(indnotuse)=[];
leftorright(indnotuse)=[];
datainfo(indnotuse,:)=[];
datanms(indnotuse)=[];

save DATA data inc azi losazienu flag_if_2D leftorright coor datainfo dem mask fault windowsize fsmpara flag_interWeight flag_adpws flag_smad  -v7.3

%% 普通加权最小二乘结果，the classical weighted least square method
Result_wls=WLS3D('DATA.mat');

dl=[-0.5,0.5];
% 画图显示三维形变，show 3-D displacements
lgdstr={'(a)E-W','(b)N-S','(c)Vertical'};
[axs,hcs]=getfig(Result_wls.enu,dl,1,1,0,0,lgdstr);
set(gcf,'Position',[133 399 2341 804],'name','3D_WLS');
reshapefigure(coor);
for i=1:length(hcs)
    clbtitle(hcs(i),'[m]',14);
    axes(axs(i));hold on
coortick(coor,[0.2,0.2]);
    if length(fault(:))~=1
        for j=1:length(fault)
            subi=lonlat2sub(fault{j},coor);
            hold on
            plot(subi(:,2),subi(:,1),'-r','linewidth',2);
        end
    end
end
nm='3D_WLS.png';
print(gcf,nm,'-dpng','-r300');
rmfigbg(nm);
pause(3);

%% SM-VCE求解三维形变， Calculate the 3-D displacements by SM-VCE
Result_smvce=SMVCE_solve3D('DATA.mat');

% Result_smvce=interpNS(Result_smvce);
save(['Result_SMVCE_',datestr(clock,'yyyymmddHHMMSS')],'Result_smvce','-v7.3');
%
%% Show Result, 画图
dl=[-0.2,0.2]*2.5;
lgdstr={'(a)E-W','(b)N-S','(c)Vertical'};
[axs,hcs]=getfig(Result_smvce.enu,dl,1,1,0,0,lgdstr);
set(gcf,'Position',[133 589 2341 614],'name','3D_SMVCE');
reshapefigure(coor);
for i=1:length(hcs)
    clbtitle(hcs(i),'[m]',14);
    axes(axs(i));hold on
coortick(coor,[0.2,0.2]);
    if length(fault(:))~=1
        for j=1:length(fault)
            subi=lonlat2sub(fault{j},coor);
            hold on
            plot(subi(:,2),subi(:,1),'-r','linewidth',2);
        end
    end
end

[row,col,~]=size(Result_smvce.enu);
nv=15;
ii=round(linspace(1,row,nv));
jj=round(linspace(1,col,nv));
[jjj,iii]=meshgrid(jj,ii);
en=Result_smvce.enu(ii,jj,1:2);
xyen=[jjj(:),iii(:),reshape(en,[],2)];
plot2d_vector(xyen,250);
% getfig(dataout);

nm='3D_SMVCE.png';
print(gcf,nm,'-dpng','-r300');
rmfigbg(nm);

ii=round(0.8*row):row;
jj=round(0.8*col):col;
enustd=Result_smvce.enu(ii,jj,:);
enustd=reshape(enustd,[],3);
enustdv=nanstd(enustd,1);

recsub=[ii(1),jj(1);
    ii(1),jj(end);
    ii(end),jj(end);
    ii(end),jj(1);
    ii(1),jj(1)];
recll=sub2lonlat(recsub,coor);

%% Show Result, 画图 方差
dl=[-3,1];
lgdstr={'(a)E-W','(b)N-S','(c)Vertical'};
varenu=log10(Result_smvce.var.enu);
[axs,hcs]=getfig(varenu,dl,1,0,0,0,lgdstr);
for i=1:length(hcs)
    clbtitle(hcs(i),'log_{10}(var)',14);
    axes(axs(i));hold on
    if length(fault(:))~=1
        for j=1:length(fault)
            subi=lonlat2sub(fault{j},coor);
            hold on
            plot(subi(:,2),subi(:,1),'-r','linewidth',2);
        end
    end
end
set(gcf,'Position',[84 281 1761 462],'name','3D_SMVCE');


nm='3D_var_SMVCE.png';
print(gcf,nm,'-dpng','-r300');
rmfigbg(nm);