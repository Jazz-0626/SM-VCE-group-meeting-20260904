function [dataout,dataout2,splitlines,ind] = quadtree_ljh(datain,var_min,num_min,ind) 
%quadtree_ljh
%--by LiuJH 2019/03/18
%--This code is used to divide images by the threshod of variance and
%window size
%
%Input:
%  datain :the input data 
%  var_min:the variance threshod of data in the final small window
%  num_min:the minmum of pixel num in the final small window
%  ind    :the global variable indicating how many small windows in the
%  last
%
%Output:
%  dataout   :with the same size of datain, but same values in each final
%  small windows
%  dataout2  :with the same size of datain, but split lines with value of
%  1, and others 0
%  splitlines:each row indicating a split line, [row1,col1,row2,col2]
%  ind       :the global variable indicating how many small windows in the
%  last
%
ind=ind+1;
dataout=datain;
dataout2=zeros(size(dataout));
splitlines=[];
if ~isempty(datain) 
    [row,col]=size(datain); 
    col2=round(col/2); 
    row2=round(row/2); 
    child1=datain(1:row2,1:col2); 
    child2=datain(1:row2,col2+1:col); 
    child3=datain(row2+1:row,1:col2); 
    child4=datain(row2+1:row,col2+1:col); 
    v=nanstd(datain(:));
    %variance > var_min and window size > num_min
    %then further divede
    if (v>var_min)&&(row*col>num_min) 
        [child1,dataout21,splitlines1,ind]=quadtree_ljh(child1,var_min,num_min,ind); 
        [child2,dataout22,splitlines2,ind]=quadtree_ljh(child2,var_min,num_min,ind); 
        [child3,dataout23,splitlines3,ind]=quadtree_ljh(child3,var_min,num_min,ind); 
        [child4,dataout24,splitlines4,ind]=quadtree_ljh(child4,var_min,num_min,ind); 
        dataout=[child1,child2;child3,child4];
        dataout2=[dataout21,dataout22;dataout23,dataout24];
        dataout2(:,col2)=1;
        dataout2(row2,:)=1;
        
        if ~isempty(splitlines2)
            splitlines2(:,2:2:4)=splitlines2(:,2:2:4)+col2;
        end
        if ~isempty(splitlines3)
            splitlines3(:,1:2:3)=splitlines3(:,1:2:3)+row2;
        end
        if ~isempty(splitlines4)
            splitlines4(:,1:2:3)=splitlines4(:,1:2:3)+row2;
            splitlines4(:,2:2:4)=splitlines4(:,2:2:4)+col2;
        end
        
        splitlines=[0,col2,row,col2;...
                    row2,0,row2,col;...
                    splitlines1;...
                    splitlines2;...
                    splitlines3;...
                    splitlines4];

    else
        if (nanmean(dataout(:))~=1)&&(~isempty(dataout))
            dataout=ones(row,col)*ind; 
            dataout(isnan(datain))=nan;
        end
    end
end
end