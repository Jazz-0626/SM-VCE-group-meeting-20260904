function [data,inc,azi,losazienu,leftorright,coor,dem,mask,fault,datainfo]=SMVCE_readdispdata()
%SMVCE_readdispdata
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to read the prepared data or the following procedure

pdata=[pwd,filesep,'SMVCE_DATA'];
pdata_information=[pdata,filesep,'data_information'];
fid=fopen(pdata_information,'r');
N=0;
datainfo={};
while ~feof(fid)
    line=fgetl(fid);
    if strcmp('#',line(1))
        continue;
    end
    N=N+1;
    datainfoi=splitstr(line);
    datainfo=[datainfo;datainfoi];
end
fclose(fid);
inc=[];
azi=[];
losazienu=zeros(N,1);
leftorright=zeros(N,1);
data=[];
for i=1:N
    losazienu(i)=str2double(datainfo{i,2});
    leftorright(i)=str2double(datainfo{i,5});
    ptifi=[pdata,filesep,datainfo{i,1}];
    datai=geotiffread(ptifi);
    data=cat(3,data,datai);
    if i==1
        datainfoi=geotiffinfo(ptifi);
        coor.corner_lon=datainfoi.RefMatrix(3,1);
        coor.corner_lat=datainfoi.RefMatrix(3,2);
        coor.post_lon=datainfoi.RefMatrix(2,1);
        coor.post_lat=datainfoi.RefMatrix(1,2);
        coor.nlines=size(datai,1);
        coor.width=size(datai,2);
    end
    [row,col]=size(datai);
    if isnan(str2double(datainfo{i,3}))
        inci=geotiffread([pdata,filesep,datainfo{i,3}]);
    else
        inci=str2double(datainfo{i,3})*ones(row,col);
    end
    inc=cat(3,inc,inci);
    if isnan(str2double(datainfo{i,4}))
        azii=geotiffread([pdata,filesep,datainfo{i,4}]);
    else
        azii=str2double(datainfo{i,4})*ones(row,col);
    end
    azi=cat(3,azi,azii);
end
[row,col,~]=size(data);
if exist([pdata,filesep,'dem.tif'],'file')
    dem=geotiffread([pdata,filesep,'dem.tif']);
else
    dem=randn(row,col);
end
if exist([pdata,filesep,'mask.tif'],'file')
    mask=geotiffread([pdata,filesep,'mask.tif']);
else
    mask=ones(row,col);
%     data1=data;
%     data1(isnan(data1))=0;
%     data1s=sum(data1~=0,3)==size(data,3);
%     [cs,rs]=meshgrid(1:col,1:row);
%     x=cs(data1s==1);
%     y=rs(data1s==1);
%     k=boundary(x,y);
%     [in,on]=inpolygon(cs,rs,x(k),y(k));
%     mask=in+on;
end
if exist([pdata,filesep,'fault.xy'],'file')
    fault=readfault([pdata,filesep,'fault.xy']);
else
    fault=0;
end
end
%% other functions

function fault=readfault(pfault)
%read exist faults data

fid=fopen(pfault,'r');
faultidex=0;

while ~feof(fid)
    line=fgetl(fid);
    if isempty(line)||strcmp(line,' ')
        continue;
    end
    if ~isempty(strfind(line,'>'))
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
