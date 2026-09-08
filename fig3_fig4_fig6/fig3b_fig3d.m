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


%%% fig3(b)
for i=1:4 % 6
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
legend('Select the telegraph model', 'Select the cross-talk pathway model', 'Select the three-state model', 'Select cross-talk three-state model')%, 'Select the positive model', 'Select negative model')

HD_F1 = [data_F1(:, 15), data_F1(:, 29)  ]; 
HD_F2 = [data_F2(:, 15), data_F2(:, 29)  ]; 
HD_E1 = [data_E1(:, 15), data_E1(:, 29)  ]; 
HD_E2 = [data_E2(:, 15), data_E2(:, 29)  ]; 
HD_F1_tw = HD_F1(AICc_F1(:,1)==AICc_F11, :);  HD_F1_c = HD_F1(AICc_F1(:,2)==AICc_F11, :);
HD_F2_tw = HD_F2(AICc_F2(:,1)==AICc_F21, :);  HD_F2_c = HD_F2(AICc_F2(:,2)==AICc_F21, :);
HD_E1_tw = HD_E1(AICc_E1(:,1)==AICc_E11, :);  HD_E1_c = HD_E1(AICc_E1(:,2)==AICc_E11, :);
HD_E2_tw = HD_E2(AICc_E2(:,1)==AICc_E21, :);  HD_E2_c = HD_E2(AICc_E2(:,2)==AICc_E21, :);
group1=[repmat({'F1_tw'}, size(HD_F1_tw(:,1),1),1); repmat({'F1_c'}, size(HD_F1_c(:,1),1),1); repmat({'F2_tw'}, size(HD_F2_tw(:,1),1),1); repmat({'F2_c'}, size(HD_F2_c(:,1),1),1); repmat({'E1_tw'}, size(HD_E1_tw(:,1),1),1); repmat({'E1_c'}, size(HD_E1_c(:,1),1),1); repmat({'E2_tw'}, size(HD_E2_tw(:,1),1),1); repmat({'E2_c'}, size(HD_E2_c(:,1),1),1) ];

%%% fig3(d)
figure(2)
boxplot([ HD_F1_tw(:,1)./HD_F1_tw(:,2); HD_F1_c(:,1)./HD_F1_c(:,2); HD_F2_tw(:,1)./HD_F2_tw(:,2); HD_F2_c(:,1)./HD_F2_c(:,2); HD_E1_tw(:,1)./HD_E1_tw(:,2); HD_E1_c(:,1)./HD_E1_c(:,2); HD_E2_tw(:,1)./HD_E2_tw(:,2); HD_E2_c(:,1)./HD_E2_c(:,2)  ], group1)
ylim([0.88, 1.56])



