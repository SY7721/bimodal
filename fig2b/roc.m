clc
clear all

kk = 1; % 1, 2, 3, 4
if kk == 1 % tw_th
     data_table = readtable( "threestate_fitmethod_bimodal_1000_1.csv"); 
     data = table2array(data_table); 
     data(:,1)=[];
     clear data_table
     data_table = readtable( "twostate_fitmethod_bimodal_1000.csv"); 
     data1 = table2array(data_table); 
     data1(:,1)=[];
elseif kk == 2 % tw_c
     data_table = readtable( "crosstalk_fitmethod_bimodal_1000.csv"); 
     data = table2array(data_table); 
     data(:,1)=[];
     clear data_table
     data_table = readtable( "twostate_fitmethod_bimodal_1000.csv"); 
     data1 = table2array(data_table); 
     data1(:,1)=[];
elseif kk == 3 % th_cth
     data_table = readtable( "crosstalk_threestate_fitmethod_bimodal_1000.csv"); 
     data = table2array(data_table); 
     data(:,1)=[];
     clear data_table
     data_table = readtable( "threestate_fitmethod_bimodal_1000_1.csv"); 
     data1 = table2array(data_table); 
     data1(:,1)=[];
elseif kk == 4 % c_cth
     data_table = readtable( "crosstalk_threestate_fitmethod_bimodal_1000.csv"); 
     data = table2array(data_table); 
     data(:,1)=[];
     clear data_table
     data_table = readtable( "crosstalk_fitmethod_bimodal_1000.csv"); 
     data1 = table2array(data_table); 
     data1(:,1)=[];
end

N1=10000;
GN=1000;
expression = cell(GN+1, N1+1);  
expression1 = cell(GN+1, N1+1);  
for i=1:GN
    expression{i+1,1}=sprintf('gene%d', i);
    expression1{i+1,1}=sprintf('gene%d', i);
end
for i = 1:N1
    expression{1, i+1} = sprintf('cell%d', i);
    expression1{1, i+1} = sprintf('cell%d', i);
end
expression(2:end, 2:end) = num2cell(NaN(GN, N1));
expression1(2:end, 2:end) = num2cell(NaN(GN, N1));

for zushu=1:1:GN
    zushu
tic
%%-------------------------------------
clearvars -except expression expression1 zushu1 data zushu filename N1 data1 filename1  kk  GN 

global  xdata tdata tt   

if zushu<=500
    if kk == 1 || kk ==2
       parameter = data1(zushu, 1:3);
       SS = telegraph(parameter, N1);
    elseif kk == 3
       parameter = data1(zushu, 1:4);
       SS = threestate(parameter, N1);
    elseif kk == 4
       parameter = data1(zushu, 1:5);
       SS = crosstalk(parameter, N1);
    end
  
elseif zushu>500
    if kk == 1
       parameter = data(zushu-500, 1:4);
       SS = threestate(parameter, N1);
    elseif kk == 2
       parameter = data(zushu-500, 1:5);
       SS = crosstalk(parameter, N1);
    elseif kk == 3 || kk ==4
       parameter = data(zushu-500, 1:6);
       SS = crosstalk_threestate(parameter, N1);
    end

end
expression(zushu+1,2:end)=num2cell(SS);

for i = 1:max(SS)+1
    y(i) = sum(SS==(i-1))./N1;
end
xdata=[]; tdata=[]; tt=[];
xdata=y;
tt=N1.*xdata;
S=numel(xdata);
for i=1:1:S;
tdata(i)=i-1;
end
clear parameter 
    if kk == 1 || kk ==2
        model = 1; 
    elseif kk == 3
        model = 3;
    elseif kk ==4
        model=2;  
    end
    tic
    parameter=MLE(model);
    toc
    if kk == 1 || kk == 2
        SS1 = telegraph(parameter(1:end-1), N1);
    elseif kk == 3
        SS1 = threestate(parameter(1:end-1), N1);
    elseif kk == 4
        SS1 = crosstalk(parameter(1:end-1), N1);
    end
expression1(zushu+1,2:end)=num2cell(SS1);

end
if kk == 1
    writecell(expression, 'tw_th_expression.csv');
    writecell(expression1, 'tw_th_expression1.csv');
elseif kk == 2
    writecell(expression, 'tw_c_expression.csv');
    writecell(expression1, 'tw_c_expression1.csv');
elseif kk == 3
    xwritecell(expression, 'th_cth_expression.csv');
    writecell(expression1, 'th_cth_expression1.csv');
elseif kk == 4
    writecell(expression, 'c_cth_expression.csv');
    writecell(expression1, 'c_cth_expression1.csv');
end
    


