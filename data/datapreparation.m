clear

% load data
USA=xlsread('data.xlsx','USA');tcode_USA=USA(1,:);k=3;
UK=xlsread('data.xlsx','UK');tcode_UK=UK(1,:);
DEU=xlsread('data.xlsx','DEU');tcode_DEU=DEU(1,:);
FRA=xlsread('data.xlsx','FRA');tcode_FRA=FRA(1,:);
ITA=xlsread('data.xlsx','ITA');tcode_ITA=ITA(1,:);
JPN=xlsread('data.xlsx','JPN');tcode_JPN=JPN(1,:);
CAN=xlsread('data.xlsx','CAN');tcode_CAN=CAN(1,:);
USA=USA(2:end,:);UK=UK(2:end,:);DEU=DEU(2:end,:);FRA=FRA(2:end,:);
ITA=ITA(2:end,:);JPN=JPN(2:end,:);CAN=CAN(2:end,:);

% transform the data based on the transformation code, where: 0: No
% transformation; 1: 400*log(y_t/y_t-1);  4. 4*y_t; 
USA_t=zeros(188,k);
for i=1:k
    if tcode_USA(i)==0
        USA_t(:,i)=USA(:,i);
    elseif tcode_USA(i)==1
        tmpz=log(USA(:,i));
        tmpz_lag=mlag2(tmpz,1);
        USA_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_USA(i)==4
        USA_t(:,i)=4*USA(:,i);        
    end
end

UK_t=zeros(188,k);
for i=1:k
    if tcode_UK(i)==0
        UK_t(:,i)=UK(:,i);
    elseif tcode_UK(i)==1
        tmpz=log(UK(:,i));
        tmpz_lag=mlag2(tmpz,1);
        UK_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_UK(i)==4
        UK_t(:,i)=4*UK(:,i);   
    end
end

DEU_t=zeros(188,k);
for i=1:k
    if tcode_DEU(i)==0
        DEU_t(:,i)=DEU(:,i);
    elseif tcode_DEU(i)==1
        tmpz=log(DEU(:,i));
        tmpz_lag=mlag2(tmpz,1);
        DEU_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_DEU(i)==4
        DEU_t(:,i)=4*DEU(:,i);   
    end
end

FRA_t=zeros(188,k);
for i=1:k
    if tcode_FRA(i)==0
        FRA_t(:,i)=FRA(:,i);
    elseif tcode_FRA(i)==1
        tmpz=log(FRA(:,i));
        tmpz_lag=mlag2(tmpz,1);
        FRA_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_FRA(i)==4
        FRA_t(:,i)=4*FRA(:,i);   
    end
end

ITA_t=zeros(188,k);
for i=1:k
    if tcode_ITA(i)==0
        ITA_t(:,i)=ITA(:,i);
    elseif tcode_ITA(i)==1
        tmpz=log(ITA(:,i));
        tmpz_lag=mlag2(tmpz,1);
        ITA_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_ITA(i)==4
        ITA_t(:,i)=4*ITA(:,i);   
    end
end

JPN_t=zeros(188,k);
for i=1:k
    if tcode_JPN(i)==0
        JPN_t(:,i)=JPN(:,i);
    elseif tcode_JPN(i)==1
        tmpz=log(JPN(:,i));
        tmpz_lag=mlag2(tmpz,1);
        JPN_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_JPN(i)==4
        JPN_t(:,i)=4*JPN(:,i);   
    end
end

CAN_t=zeros(188,k);
for i=1:k
    if tcode_CAN(i)==0
        CAN_t(:,i)=CAN(:,i);
    elseif tcode_CAN(i)==1
        tmpz=log(CAN(:,i));
        tmpz_lag=mlag2(tmpz,1);
        CAN_t(:,i)=(tmpz-tmpz_lag)*400;
    elseif tcode_CAN(i)==4
        CAN_t(:,i)=4*CAN(:,i);   
    end
end

% Construct final dataset 
USA_s=USA_t(2:end,:);
UK_s=UK_t(2:end,:);
DEU_s=DEU_t(2:end,:);
FRA_s=FRA_t(2:end,:);
ITA_s=ITA_t(2:end,:);
JPN_s=JPN_t(2:end,:);
CAN_s=CAN_t(2:end,:);
data_panel_temp=[CAN_s DEU_s FRA_s ITA_s JPN_s UK_s USA_s];
[panel_temp,~]=remove_outliers(data_panel_temp);
[~,~,~,~,data_panel] = factors_em(panel_temp,8,2,2);


save('dataQ','data_panel')
