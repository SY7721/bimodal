clc
clear all
close all

[data,name]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");
[~,name1] = xlsread("expression of bimodal genes.xlsx", 'E c57');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
data_E1=data(row,:);  name_E1=name(row+1,1);
AICc_E1 = [data_E1(:,10) data_E1(:,24) data_E1(:,37) data_E1(:,52) ]; %  data(:,65) data(:,78)
AICc_E11 = min(AICc_E1,[],2);
clear data name row  name1 

[data,name]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");
[~,name1] = xlsread("expression of bimodal genes.xlsx", 'E cast');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
data_E2=data(row,:);  name_E2=name(row+1,1);
AICc_E2 = [data_E2(:,10) data_E2(:,24) data_E2(:,37) data_E2(:,52)];  
AICc_E21 = min(AICc_E2,[],2);
clear data name row  name1

[data,name]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");
[~,name1] = xlsread("expression of bimodal genes.xlsx", 'F c57');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
data_F1=data(row,:);  name_F1=name(row+1,1);
AICc_F1 = [data_F1(:,10) data_F1(:,24) data_F1(:,37) data_F1(:,52)];  AICc_F11 = min(AICc_F1,[],2);
clear data name row name1

[data,name]=xlsread("Fibroblasts_cast_fitmethod_nofixp.csv");
[~,name1] = xlsread("expression of bimodal genes.xlsx", 'F cast');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
data_F2=data(row,:);  name_F2=name(row+1,1);
AICc_F2 = [data_F2(:,10) data_F2(:,24) data_F2(:,37) data_F2(:,52)]; 
AICc_F21 = min(AICc_F2,[],2);
clear data name row  name1

%%%% figS4
for i=1: 6
    F1(i)=length(find(AICc_F1(:,i)==AICc_F11))./length(AICc_F11);
    F2(i)=length(find(AICc_F2(:,i)==AICc_F21))./length(AICc_F21);
    E1(i)=length(find(AICc_E1(:,i)==AICc_E11))./length(AICc_E11);
    E2(i)=length(find(AICc_E2(:,i)==AICc_E21))./length(AICc_E21);
end

F=[F1; F2; E1; E2];
x=[1, 2, 3, 4];
h=bar(x,F,1);
set(gca,'XTickLabel',{'Fibroblasts c57', 'Fibroblasts cast', 'Embryonic c57', 'Embryonic cast'}) 
ylim([0, 1])
legend('Select the telegraph model', 'Select the cross-talk pathway model', 'Select the three-state model', 'Select cross-talk three-state model', 'Select the positive model', 'Select negative model')



clc
clear all

[data,name]=xlsread("Fibroblasts_c57_AICc_bimodal4.xlsx");
[data1,name1]=xlsread("Fibroblasts_cast_AICc_bimodal4.xlsx");
[data2,name2]=xlsread("Embryonic_c57_AICc_bimodal4.xlsx");
[data3,name3]=xlsread("Embryonic_cast_AICc_bimodal4.xlsx");
n=[224, 224, 188, 188]; 
i=1; k=4;
Smin_min31=data(:,72)-2.*k-2*k*(k+1)./(n(i)-k-1);
k=5;
three1=2*k+Smin_min31+2*k*(k+1)./(n(i)-k-1);
i=2; k=4;
Smin_min32=data1(:,72)-2.*k-2*k*(k+1)./(n(i)-k-1);
k=5;
three2=2*k+Smin_min32+2*k*(k+1)./(n(i)-k-1);
i=3; k=4;
Smin_min33=data2(:,72)-2.*k-2*k*(k+1)./(n(i)-k-1);
k=5;
three3=2*k+Smin_min33+2*k*(k+1)./(n(i)-k-1);
i=4; k=4;
Smin_min34=data3(:,72)-2.*k-2*k*(k+1)./(n(i)-k-1);
k=5;
three4=2*k+Smin_min34+2*k*(k+1)./(n(i)-k-1);

AICc_F1=[data(:,23) data(:,39) data(:,55)  three1 ];
AICc_F2=[data1(:,23) data1(:,39) data1(:,55) three2];
AICc_E1=[data2(:,23) data2(:,39) data2(:,55)  three3];
AICc_E2=[data3(:,23) data3(:,39) data3(:,55)  three4];

AICc=min(AICc_F1,[],2);  AICc1=min(AICc_F2,[],2);  AICc2=min(AICc_E1,[],2);  AICc3=min(AICc_E2,[],2);

F1=[length(find(AICc==AICc_F1(:,1))) length(find(AICc==AICc_F1(:,2)))  length(find(AICc==AICc_F1(:,3))) length(find(AICc==AICc_F1(:,4))) ]./length(AICc);
F2=[length(find(AICc1==AICc_F2(:,1))) length(find(AICc1==AICc_F2(:,2)))  length(find(AICc1==AICc_F2(:,3))) length(find(AICc1==AICc_F2(:,4)))  ]./length(AICc1);
E1=[length(find(AICc2==AICc_E1(:,1))) length(find(AICc2==AICc_E1(:,2)))  length(find(AICc2==AICc_E1(:,3))) length(find(AICc2==AICc_E1(:,4)))]./length(AICc2);
E2=[length(find(AICc3==AICc_E2(:,1))) length(find(AICc3==AICc_E2(:,2)))  length(find(AICc3==AICc_E2(:,3))) length(find(AICc3==AICc_E2(:,4))) ]./length(AICc3);

F=[F1; F2; E1; E2];
x=[1, 2, 3, 4];
h=bar(x,F,1);
set(gca,'XTickLabel',{'Fibroblasts c57', 'Fibroblasts cast', 'Embryonic c57', 'Embryonic cast'})
ylim([0, 1])
legend('Select the two-state model', 'Select the cross-talk model', 'Select the three-state model', 'Select cross-talk three-state model')

