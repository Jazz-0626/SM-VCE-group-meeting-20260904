function Bgeo=lk_vec(azi,inc,losazi,leftright)
%lk_vec
%--By Liu Jihong 2021/09/26 @ Central South University
%--This function is used to calculate the coefficient between InSAR and 3-D displacements
%
%Modified:
%20211122: add the coefficient calculation for ew, ns, vertical displacement
[row,col,data_num]=size(azi);
Bgeo=zeros(row,col,3*data_num);
for i=1:data_num
    switch losazi(i)
        case 1
            Bgeo(:,:,1+3*(i-1))=-leftright(i)*sind(inc(:,:,i)).*sind(azi(:,:,i)-270);
            Bgeo(:,:,2+3*(i-1))=-leftright(i)*sind(inc(:,:,i)).*cosd(azi(:,:,i)-270);
            Bgeo(:,:,3*i)=cosd(inc(:,:,i));
        case 2
            Bgeo(:,:,1+3*(i-1))=-cosd(azi(:,:,i)-270);
            Bgeo(:,:,2+3*(i-1))=sind(azi(:,:,i)-270);
            Bgeo(:,:,3*i)=zeros(row,col);
        case 3
            Bgeo(:,:,1+3*(i-1))=1;
            Bgeo(:,:,2+3*(i-1))=0;
            Bgeo(:,:,3*i)=0;
        case 4
            Bgeo(:,:,1+3*(i-1))=0;
            Bgeo(:,:,2+3*(i-1))=1;
            Bgeo(:,:,3*i)=0;
        case 5
            Bgeo(:,:,1+3*(i-1))=0;
            Bgeo(:,:,2+3*(i-1))=0;
            Bgeo(:,:,3*i)=1;
    end
end
if row*col*data_num==1
    Bgeo=reshape(Bgeo,1,3);
end
end