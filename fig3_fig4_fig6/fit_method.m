
clc
clear all

nn = 2; % 1: F c57; 2: F cast; 3: E c57; 4: E c57
if nn == 1
    raw = readcell('Fibroblasts_c57_distribution.csv');
    filename = 'Fibroblasts_c57_fitmethod_nofixp.csv';
elseif nn == 2
    raw = readcell('Fibroblasts_cast_distribution.csv');
    filename = 'Fibroblasts_cast_fitmethod_nofixp.csv';
elseif nn == 3
    raw = readcell('Embryonic_c57_distribution.csv');
    filename = 'Embryonic_c57_fitmethod_nofixp.csv';
elseif nn == 4
    raw = readcell('Embryonic_cast_distribution.csv');
    filename = 'Embryonic_cast_fitmethod_nofixp.csv';
end
name = raw(1,:);
data = raw(2:end,:);
name2 = name;
data2 = str2double(string(data));
clear raw data name 

mmm1 = [224, 224, 188, 188]; 
CN=mmm1(nn); 
[~,GN]=size(data2);
clear mmm1 
str={'name', 'mean', 'fano', 'skewness', 'kurtosis', 'p+telegraph',  'la', 'ga', 'v', 'p',  'AICc1',  'mean1', 'fano1', 'skewness1', 'kurtosis1', 'HD1', 'type1',  'p+cross-talk','la1','la2','q1', 'gamma', 'nu', 'p', 'AICc2', 'mean2', 'fano2', 'skewness2',  'kurtosis2', 'HD2', 'type2', 'p+threestate','la1','la2','gamma', 'nu', 'p', 'AICc2', 'mean2', 'fano2', 'skewness2',  'kurtosis2', 'HD2', 'type2', 'p+three-sate cross-talk','la1','kappa1', 'kappa2', 'q1', 'gamma', 'nu', 'p', 'AICc2', 'mean2', 'fano2', 'skewness2',  'kurtosis2', 'HD2', 'type2', 'p+positive',  'la', 'ga', 'v', 'mu', 'p',  'AICc1',  'mean1', 'fano1', 'skewness1', 'kurtosis1', 'HD1', 'type1', 'p+negative',  'la', 'ga', 'v', 'nu', 'p',  'AICc1',  'mean1', 'fano1', 'skewness1', 'kurtosis1', 'HD1', 'type1'};
n_rows = GN;
value_final=cell(n_rows+1, length(str));
value_final(1, :) = str;
value_final(2:end, 2:end) = num2cell(NaN(n_rows, length(str)-1));
value_final(2:end,1) = name2(1:end);

for zushu =1: 1: GN
    %%------------------------------------
    clearvars -except zushu  filename GN CN sheet data2 name2 value_final nn
    global  xdata tt tdata 
    zushu
    
    xdata=[]; tdata=[]; tt=[];
    xdata=clearnan(data2(5:end,zushu));
    S=length(xdata);
    for i = 1:1:S;
        tdata(i) = i-1;
        tt(i) = CN.*xdata(i);
    end
    %%% moment
    mean=0; twom1=0; threem=0; fourm=0;
    for t=1:1:S
        mean=mean+(t-1).*xdata(t);
        twom1=twom1+((t-1)^2).*xdata(t);
        threem=threem+((t-1)^3).*xdata(t);
        fourm=fourm+((t-1)^4).*xdata(t);
    end
    fano=twom1/mean-mean;
    Var=0; threec=0;fourc=0;
    for i=1:1:S
        Var=Var+((i-1-mean)^2).*xdata(i);
        threec=threec+((i-1-mean)^3).*xdata(i);
        fourc=fourc+((i-1-mean)^4).*xdata(i);
    end
    skewness=threec/(Var^(3/2));     kurtosis=fourc/(Var^2)-3;
    value_final(zushu + 1,2:5) = num2cell([ mean  fano   skewness kurtosis]);
    clear mean fano  skewness kurtosis Var threec fourc twom1 threem fourm
    
    models = [1, 2, 3, 4, 5, 6]; %% 1: Binomial +telelgraph; 2: Binomial + crosstalk; 3: Binomial + threestate; 4:Binomial + crosstalk-threestate; 5: Binomial + positive; 6: Binomial + negative
    k_values = [4, 5, 5, 6, 5, 5];
    start_cols = [7, 19, 33, 46, 61, 74];
    start_cols1 = [17, 31, 44, 59, 72, 85];
for m_idx = 1:length(models)
    model = models(m_idx);
    k = k_values(m_idx);
    start_col = start_cols(m_idx);
    start_col1 = start_cols1(m_idx);
    tic
    parameter = MLE_p(model);
    toc
    Smin_min = parameter(end);
    AICc = caculate_AICc(Smin_min, k, CN);
    parameter1 = parameter(1:end-1);
    [mean_val, fano_val, skewness_val, kurtosis_val, HD, type] = caculate_moment_p(parameter1, model);
   
    value_final(zushu + 1, start_col:start_col1) =  num2cell([parameter1, AICc, mean_val, fano_val, skewness_val, kurtosis_val, HD, type]);
    
    clear model parameter parameter1 AICc mean_val fano_val skewness_val kurtosis_val HD type k Smin_min
end
   
    zushu
     
end
writecell(value_final, filename);


