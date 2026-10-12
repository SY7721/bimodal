clc
clear all

data_table = readtable( "crosstalk_threestate_fitmethod_bimodal_1000.csv"); % crosstalk threestate
data = table2array(data_table); 
data(:,1)=[];
clear data_table
data_table = readtable( "crosstalk_fitmethod_bimodal_1000.csv"); % threestate
data1 = table2array(data_table); 
data1(:,1)=[];

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

for zushu=1:1:500  
    zushu
tic
%%-------------------------------------
clearvars -except expression expression1 zushu1 data zushu filename N1 data1 filename1  kkkk  GN 

global  xdata tdata tt   

if zushu<=500
   parameter = data1(zushu, 1:5);
   SS = crosstalk(parameter, N1);
  
elseif zushu>500
     parameter = data(zushu-500, 1:6);
     SS = crosstalk_threestate(parameter, N1);

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

    model=2;  
    tic
    parameter=MLE(model);
    toc

    SS1 = crosstalk(parameter(1:end-1), N1);
expression1(zushu+1,2:end)=num2cell(SS1);

end

writecell(expression, 'c_cth_expression.csv');
writecell(expression1, 'c_cth_expression1.csv');