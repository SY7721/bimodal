
clc
clear all
close all

N1 = [100 250 500 1000 10000];
kk = 1;
CN= N1(kk); %% cell number
GN=500; 
model = 1; %% 1: two-state; 2: crosstalk; 3: threestate; 4:crosstalk-threestate
if model == 1
    data = xlsread("twostate_fitmethod_bimodal_infty.xlsx");
    filename = 'twostate_bimodal_fit.xlsx';
    filename1 = 'twostate_bimodal_expression.xlsx';
    str = {'real twostate data', 'lam', 'ga', 'v', 'mean', 'fano', 'Skewness', 'Kurtosis',  '矩方法', 'la', 'ga', 'v', 'mean', 'fano',  'Skewness', 'Kurtosis','HD','type', '两状态估计','la', ' ga','  v', ' AICc','mean', 'fano',  'Skewness', 'Kurtosis', 'HD', 'type', '两路径估计','la1', 'la2','q1', 'ga',' v', ' AICc_5',' mean','fano','Skewness', 'Kurtosis', 'HD' , 'type','三状态估计','la1', 'la2', ' ga','  v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type', 'crosstalk-threestate估计', 'lam','kappa1', 'kappa2', 'q1',  'ga', 'v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type'    };
elseif model == 2
    data = xlsread("crosstalk_fitmethod_bimodal_infty.xlsx");
    filename = 'crosstalk_bimodal_fit.xlsx';
    filename1 = 'crosstalk_bimodal_expression.xlsx';
    str = {'real cross-talk data', 'lam1','lam2', 'q1',  'ga', 'v', 'mean', 'fano', 'Skewness', 'Kurtosis',  '矩方法', 'la', 'ga', 'v', 'mean', 'fano',  'Skewness', 'Kurtosis','HD','type', '两状态估计','la', ' ga','  v', ' AICc','mean', 'fano',  'Skewness', 'Kurtosis', 'HD', 'type', '两路径估计','la1', 'la2','q1', 'ga',' v',' AICc_5',' mean','fano','Skewness', 'Kurtosis', 'HD' , 'type','三状态估计','la1', 'la2', ' ga','  v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type', 'crosstalk-threestate估计', 'lam','kappa1', 'kappa2', 'q1',  'ga', 'v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type'    };
elseif model == 3
    data = xlsread("threestate_fitmethod_bimodal_infty.xlsx");
    filename = 'threestate_bimodal_fit.xlsx';
    filename1 = 'threestate_bimodal_expression.xlsx';
    str = {'real three-state data', 'lam1','lam2', 'ga', 'v', 'mean', 'fano', 'Skewness', 'Kurtosis',  '矩方法', 'la', 'ga', 'v', 'mean', 'fano',  'Skewness', 'Kurtosis','HD','type', '两状态估计','la', ' ga','  v', ' AICc','mean', 'fano',  'Skewness', 'Kurtosis', 'HD', 'type', '两路径估计','la1', 'la2','q1', 'ga',' v', ' AICc_5',' mean','fano','Skewness', 'Kurtosis', 'HD' , 'type','三状态估计','la1', 'la2', ' ga','  v',  'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type', 'crosstalk-threestate估计', 'lam','kappa1', 'kappa2', 'q1',  'ga', 'v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type'    };
elseif model == 4
    data = xlsread("crosstalk-threestate_fitmethod_bimodal_infty.xlsx");
    filename = 'crosstalk_threestate_bimodal_fit.xlsx';
    filename1 = 'crosstalk_threestate_bimodal_expression.xlsx';
    str = {'real cross-talk three-state data', 'lam1','kappa1','kappa2', 'q1',  'ga', 'v', 'mean', 'fano', 'Skewness', 'Kurtosis',  '矩方法', 'la', 'ga', 'v', 'mean', 'fano',  'Skewness', 'Kurtosis','HD','type', '两状态估计','la', ' ga','  v', ' AICc','mean', 'fano',  'Skewness', 'Kurtosis', 'HD', 'type', '两路径估计','la1', 'la2','q1', 'ga',' v', ' AICc_5',' mean','fano','Skewness', 'Kurtosis', 'HD' , 'type','三状态估计','la1', 'la2', ' ga','  v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type', 'crosstalk-threestate估计', 'lam','kappa1', 'kappa2', 'q1',  'ga', 'v', 'AICc',  'mean', 'fano',  'Skewness', 'Kurtosis','HD' , 'type'    };
end
NN=length(str);
value_final = cell(GN+1, NN);
expression = cell(GN+1, CN+1);  
value_final(1, :) = str;
for i=1:GN
     value_final{i+1,1}=sprintf('gene%d', i);
    expression{i+1,1}=sprintf('gene%d', i);
end
 value_final(2:end, 2:end) = num2cell(NaN(GN, NN-1));
for i = 1:CN
    expression{1, i+1} = sprintf('cell%d', i);
end
expression(2:end, 2:end) = num2cell(NaN(GN, CN));

for zushu =1 : GN
    clearvars -except zushu kk data CN N1 filename1 value_final expression GN data NN model filename
    global  xdata tt tdata 
    zushu

    if model==1
    % lam=0.1+4.9*rand(); gamma=0.1+4.9*rand(); nu=5+45*rand(); 
    lam=data(zushu,1); gamma=data(zushu,2); nu=data(zushu,3);
     parameter = [lam gamma nu]; 
     ydata = SSA_model(parameter, model, CN);
    elseif model == 2
     % lam1=0.1; lam2=0.1+4.9*rand(); q1 = rand(); gamma=0.1+4.9*rand(); nu=5+45*rand(); 
     lam1=data(zushu,1); lam2=data(zushu,2); q1 = data(zushu,3); gamma=data(zushu,4); nu=data(zushu,5); 
     parameter = [lam1 lam2 q1 gamma nu]; 
     ydata = SSA_model(parameter, model, CN);
    elseif model == 3
     %lam1=0.1+4.9*rand(); lam2=lam1; gamma=0.1+4.9*rand(); nu=5+45*rand(); 
     lam1=data(zushu,1); lam2=data(zushu,2); gamma=data(zushu,3); nu=data(zushu,4);
     parameter = [lam1 lam2 gamma nu]; 
     ydata = SSA_model(parameter, model, CN);
    elseif model == 4
     % lam=0.1+4.9*rand(); kappa1 = 0.2; kappa2 = 0.2; q1= rand(); gamma=0.1+4.9*rand(); nu=5+45*rand(); 
     lam=data(zushu,1); kappa1 = data(zushu,2); kappa2 =data(zushu,3); q1= data(zushu,4); gamma=data(zushu,5); nu=data(zushu,6); 
     parameter = [lam kappa1 kappa2 q1 gamma nu]; 
     ydata = SSA_model(parameter, model, CN);
    end

     expression(zushu+1,2:end)=num2cell(ydata);

     M = 3*ceil(nu)+1;
     for i = 1:M
         yyy(i) = sum(ydata==(i-1))/CN;
     end
     xdata=[]; tdata=[]; tt=[];
     xdata=yyy;
     tt=CN.*xdata;
     S=numel(xdata);
     for i=1:1:S;
         tdata(i)=i-1;
     end
    %%% moment
    mean_val=0; twom1=0; threem=0; fourm=0;
    for t=1:1:S
        mean_val=mean_val+(t-1).*xdata(t);
        twom1=twom1+((t-1)^2).*xdata(t);
        threem=threem+((t-1)^3).*xdata(t);
        fourm=fourm+((t-1)^4).*xdata(t);
    end
    fano=twom1/mean_val-mean_val;
    Var=0; threec=0;fourc=0;
    for i=1:1:S
        Var=Var+((i-1-mean_val)^2).*xdata(i);
        threec=threec+((i-1-mean_val)^3).*xdata(i);
        fourc=fourc+((i-1-mean_val)^4).*xdata(i);
    end
    skewness=threec/(Var^(3/2));     kurtosis=fourc/(Var^2)-3;

    end_col1 = [8, 10, 9, 11];
    value_final(zushu+1,2:end_col1(model))=num2cell([parameter mean_val fano skewness kurtosis]);
clear parameter
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 矩方法

     e1=0;e2=0;e3=0;
for i=1:1:S;
    e1=(i-1).*xdata(i)+e1;
    e2=((i-1).^2).*xdata(i)-(i-1).*xdata(i)+e2;
    e3=((i-1).^3).*xdata(i)-3*((i-1).^2).*xdata(i)+2*(i-1).*xdata(i)+e3;
end
    e1;  %一阶矩
    e2;  %二阶矩
    e3;  %三阶矩
    r1=e1;r2=e2/e1;r3=e3/e2;
    lain1=2*r1*(r3-r2)/(r1*r2-2*r1*r3+r2*r3);
    gain1=2*(r2-r1)*(r1-r3)*(r3-r2)/((r1*r2-2*r1*r3+r2*r3)*(r1-2*r2+r3));
    vin11=(-r1*r2+2*r1*r3-r2*r3)/(r1-2*r2+r3);

    value_final(zushu+1, end_col1(model)+2 : end_col1(model)+4)=num2cell([lain1 gain1 vin11]);

    if lain1>0 && gain1>0 && vin11>0
       parameter=[lain1 gain1 vin11];
       [mean5,  fano5,  skewness5, kurtosis5,HD5, type5]=caculate_moment(parameter, 1);   
       
       value_final(zushu+1, end_col1(model)+5 : end_col1(model)+10)=num2cell([ mean5 fano5 skewness5 kurtosis5 HD5 type5]);
    end
    clear parameter
    
    model1 = [1 2 3 4];
    k = [3, 4, 4, 5];
    start_col = [end_col1(model)+12 end_col1(model)+23 end_col1(model)+36 end_col1(model)+48];
    end_col = [end_col1(model)+21 end_col1(model)+34 end_col1(model)+46 end_col1(model)+60];
    for alpha = 1:length(model1)
         m = model1(alpha);  
         tic
         parameter = MLE(m);
         toc
         Smin_min=parameter(end);
         AICc = caculate_AICc(Smin_min,k(alpha),CN);
         [mean_val, fano,  Skewness, Kurtosis, HD, type] = caculate_moment(parameter(1:end-1),m);
         value_final(zushu+1, start_col(alpha) : end_col(alpha))=num2cell([parameter(1:end-1) AICc mean_val fano skewness kurtosis HD type]);
         clear parameter 
    end

end

sheet = {'100', '250', '500', '1000', '1000'};
writecell(value_final, filename, 'Sheet', sheet{kk});
writecell(expression, filename1, 'Sheet', sheet{kk});




