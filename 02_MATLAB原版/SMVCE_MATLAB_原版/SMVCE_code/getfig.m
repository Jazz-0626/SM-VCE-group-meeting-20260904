function [axs,hcs]=getfig(data,dls,flag_clb,flag_tklb,siz,newind,lgdstr)
%getfig
%--By Liu Jihong, Hu Jun, Li Zhiwei 2021/09/26 @ Central South University
%--By 刘计洪, 胡俊, 李志伟 2021/09/26 @ 中南大学
%
%--This function is used to show figures with more beautiful style, without
%much more lags
%
%Inputs:
%  data: a matrix with size of m*n*p,echo floor represents a figure
%  dls : the data limits of each figure, with size of p*2
%  flag_clb: if plot the colorbar, 1 yes (default) /0 no
%  flag_tklb: if show the ticklabels, 1 yes (default) /0 no
%  siz: the subplot figures size
%  newind: the new order of the figures
%  lgdstr：default is a,b,c,..., you can set the new ones with a string cell
%
%Example:
%  data=repmat(peaks(100),1,1,9);
%  dls=repmat([-5,5],9,1);
%  axs=getfig(data,dls);
%  axes(axs(2))
%  hold on
%  plot(10,10,'^m','markerfacecolor','m');

[row,col,figNum]=size(data);
r=floor(sqrt(figNum));
c=ceil(figNum/r);

%the data limits
if nargin==1||length(dls)==1
    dls=zeros(figNum,2);
    for i=1:figNum
        datai=data(:,:,i);
        %         datai(abs(datai)==inf)=nan;
        dls(i,:)=[nanmin(datai(:)),nanmax(datai(:))];
        a=dls(i,:);
        if sum(isnan(a))>0
            dls(i,:)=[-1,1];
        end
        if dls(i,1)==dls(i,2)
            dls(i,:)=[dls(i,1)-0.5,dls(i,1)+0.5];
        end
        if isnan(dls(i,2))||abs(sum(dls(i,:)))==inf
            dls(i,:)=[-0.5,0.5];
        end
    end
end
if nargin<=2
    flag_clb=0;
    flag_tklb=0;
end
if nargin<=3
    flag_tklb=0;
end
if nargin>=5&&length(siz)==2
    r=siz(1);
    c=siz(2);
end
if nargin<6||length(newind)==1
    newind=1:figNum;
end
if size(dls,1)==1&&figNum>1
    dls=repmat(dls,figNum,1);
end
for i=1:figNum
    if dls(i,1)==dls(i,2)&&dls(i,1)==0
        datai=data(:,:,i);
        %         datai(abs(datai)==inf)=nan;
        dls(i,:)=[nanmin(datai(:)),nanmax(datai(:))];
        if dls(i,1)==dls(i,2)
            dls(i,:)=[dls(i,1)-0.5,dls(i,1)+0.5];
        end
        if isnan(dls(i,2))||abs(sum(dls(i,:)))==inf
            dls(i,:)=[-0.5,0.5];
        end
    end
end

cbw=6/77.5/c;%the width of colorbar
tklbw=6/77.5/c;%the width of ticklabel width
cbtexw=5/77.5/c;%the width of colorbar tex
cblinter=2.5/77.5/c;%the interval width left of colorbar
lsinter=2/77.5/c;%the leftest start interval
lowerestinter=1.5/57.5/r;%the lowerest start interval
tklbh=3/57.5/r;%the width of ticklabel hight


scnsize = get(0,'MonitorPosition');%the size of screen
scnsize=scnsize(1,:);

tw=1000; %the whole image width
twp=1-tklbw*flag_tklb-cbtexw*flag_clb-flag_clb*cblinter-flag_clb*cbtexw-lsinter*2;%rm the clb width effect
th=round(tw/c*row/col*r*twp); %the whole image hight
if th>=scnsize(4)*0.9
    tw=tw*scnsize(4)*0.9/th; %the whole image width
    twp=1-tklbw*flag_tklb-cbtexw*flag_clb-flag_clb*cblinter-flag_clb*cbtexw-lsinter*2;%rm the clb width effect
    th=scnsize(4)*0.9; %the whole image hight
end
    
figw=tw/c; %each figure width
figh=round(0.9*figw*row/col);
abcfs=20*(figw+figh)/2/300;%abcfontsize
tkfs=10*(figw+figh)/2/300;%ticklabel fontsize
% if abcfs<12
%     abcfs=12;
% end
if tkfs<12
%     tkfs=12;
end

fig1=figure;
set(fig1,'Position',[scnsize(3)/2-tw/2,scnsize(4)/2-th/2,tw,th],'papertype','a1');%let the figure in the center of screen

iw=((1-lsinter*2)/c-flag_clb*(cblinter+cbw+cbtexw)-flag_tklb*(tklbw))*0.95;%the width of each image
ih=((1-lowerestinter*2)/r-flag_tklb*tklbh)*0.9;%the height
axs=[];
hcs=[];
k=0;
newind1=sort(newind,'ascend');

for i=newind
    k=k+1;
    [lie,hang]=ind2sub([c,r],i);
    
    l=(lie-1)/c*(1-lsinter)+lsinter;
    b=(r-hang)/r*(1-lowerestinter)+lowerestinter;
    if flag_tklb==1
        l=(lie-1)/c+lsinter+tklbw;
        b=(r-hang)/r+lowerestinter+tklbh;
    end
    axs(k)=axes('Position',[l b iw ih]);
        a=dls(k,:);
        if sum(isnan(a))>0
            dls(k,:)=[-1,1];
        end
    if dls(k,1)>dls(k,2)
        dls(k,1)=dls(k,2)+1;
    end
    hi=imagesc(data(:,:,k),dls(k,:));
    set(hi,'alphadata',~isnan(data(:,:,k)));
    str=sprintf('(%s)',char(find(newind1==i)+96));
    if exist('lgdstr','var')
        str=lgdstr{k};
        str=strrep(str,'_','\_');
    end
    fz=get(gca,'FontSize');
    if figNum>1
        text(0.56*15/col,0.56*15/col,str,'BackgroundColor',[1 1 1],'fontweight','bold','FontSize',fz*1.5,'horizontalAlignment','left','verticalAlignment','top');
    end
    if flag_clb==1
        hc=colorbar('Position',[l+iw+cblinter,b,cbw,ih],'FontSize',fz*1.2);%colorbar();
        set(hc,'fontsize',fz*1.2);
        hcs=[hcs,hc];
    end
    set(gca,'FontSize',fz*1.2);
    if flag_tklb==0
        set(gca,'xticklabel',[],'yticklabel',[]);
    end
end
end
