clc
clear all
close all

%%% fig6b
nn = 1; %% 1: F c57; 2: F cast; 3: E c57; 4: E cast
if nn == 1
    [data,name]=xlsread("Fibroblasts_c57_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'F c57');
    nname = 'F c57';
elseif nn == 2
    [data,name]=xlsread("Fibroblasts_cast_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'F cast');
    nname = 'F cast';
elseif nn == 3
    [data,name]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'E c57');
    nname = 'E c57';
elseif nn == 4
    [data,name]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");
    [~,name1] = xlsread("expression of bimodal genes.xlsx", 'E cast');
    nname = 'E cast';
end
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
data_F1=data(row,:);  name_F1=name(row+1,1);
AICc_F1 = [data_F1(:,10) data_F1(:,24) ]; %  data(:,65) data(:,78)
AICc_F11 = min(AICc_F1,[],2);
type=[data_F1(:,16) data_F1(:,30) ]; 
data_tw = data_F1(AICc_F1(:,1) == AICc_F11, :); name_tw = name_F1(AICc_F1(:,1) == AICc_F11, :); 
data_c =  data_F1(AICc_F1(:,2) == AICc_F11, :);  name_c = name_F1(AICc_F1(:,2) == AICc_F11, :);  

[data,name]=xlsread("Embryonic_cast_distribution.csv");
[coom1,~, idx1] = intersect(name_tw, name, 'stable');
data1 = data(:, idx1);
[coom2, ~, idx2] = intersect(name_c, name, 'stable');
data2 = data(:, idx2);

SS = length(data1); kb1 = nan(SS,1);
for zushu=1: SS  
    xdata=[];
    xdata = clearnan(data1(5:end, alpha));
  
   S=length(xdata);
   parameter=data_tw(zushu,6:8);
   M=S-1;
   model = 1;
   dist = fsp_model(parameter, model, M);
   x=dist;

 hhigh =NaN; hlow =NaN; valley =NaN;
 hhigh=x(1);
for i=3:N-1
    if x(i-1)<x(i) && x(i)>=x(i+1)       
        hlow = min([hhigh,x(i)]);
        hhigh = max([hhigh,x(i)]);
    end 
end
for i=2:N-1
    if x(i)<x(i-1) && x(i)<x(i+1)
        valley = x(i);
    end
end
if all(~isnan([hlow valley hhigh]))
   kb1(zushu)=(hlow-valley)/hhigh;
end
 clear xdata S M N parameter x dist Qmod ssd vec num one temp Q 
end

clear data2 name2 is_member data3 SS 

SS = length(data2); kb2 = nan(SS,1);
for zushu=1: SS  
   
   xdata=[];
   xdata = clearnan(data2(5:end, alpha));  
   S=length(xdata);
   parameter=data_c(zushu,18:22);
   M=S-1;
   model = 2;
   dist = fsp_model(parameter, model, M);
   x=dist;
  hhigh =NaN; hlow =NaN; valley =NaN;
 hhigh=x(1);
for i=3:N-1
    if x(i-1)<x(i) && x(i)>=x(i+1)       
        hlow = min([hhigh,x(i)]);
        hhigh = max([hhigh,x(i)]);
    end 
end
for i=2:N-1
    if x(i)<x(i-1) && x(i)<x(i+1)
        valley = x(i);
    end
end
if all(~isnan([hlow valley hhigh]))
   kb2(zushu)=(hlow-valley)/hhigh;
end
clear xdata S M N parameter x dist Qmod ssd vec num one temp Q 
end

group1 = [
    repmat({'twostate'}, size(clearnan(kb1),1), 1);   
    repmat({'crosstalk'}, size(clearnan(kb2),1), 1)];
figure(1)
boxplot([clearnan(kb1); clearnan(kb2)], group1)
ylim([0 0.52])
ylabel('bimodal strength')
title('nname')


