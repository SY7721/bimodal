clc
clear all

sample = 1; % 1: 100; 2: 250; 3: 500; 4: 1000; 5:10000
sheet = {'100', '250', '500', '1000', '10000'};
[data,name]=xlsread("twostate_bimodal_fit.xlsx", sheet{sample});
[data1,name1]=xlsread("threestate_bimodal_fit.xlsx", sheet{sample});
[data2,name2]=xlsread("crosstalk_bimodal_fit.xlsx", sheet{sample});
[data3,name3]=xlsread("crosstalk_threestate_bimodal_fit.xlsx", sheet{sample});

%%%% fig2(d)
kurt_r2=data(:,7); kurt_m2=data(:,15);
kurt_r3=data1(:,8); kurt_m3=data1(:,16);
kurt_rc=data2(:,9); kurt_mc=data2(:,17);
kurt_r2_u=data3(:,10); kurt_m2_u=data3(:,18);

kurt2=clearnan(kurt_r2-kurt_m2);
kurt3=clearnan(kurt_r3-kurt_m3);
kurt_c=clearnan(kurt_rc-kurt_mc);
kurt2_u=clearnan(kurt_r2_u-kurt_m2_u);

group1=[repmat({'twostate'}, size(kurt2,1),1);  repmat({'threestate'}, size(kurt3,1),1);  repmat({'crosstalk'}, size(kurt_c,1),1);  repmat({'cross-talk three-state'}, size(kurt2_u,1),1)];
figure(1)
boxplot([kurt2;  kurt3;  kurt_c;  kurt2_u], group1)
title(sheet{sample})
clear group1

%%%% fig2(c)
if sample == 5
    [data01,~]=xlsread("twostate_bimodal_expression.xlsx", sheet{sample});
    [data11,~]=xlsread("threestate_bimodal_expression.xlsx", sheet{sample});
    [data21,~]=xlsread("crosstalk_bimodal_expression.xlsx", sheet{sample});
    [data31,~]=xlsread("crosstalk_threestate_bimodal_expression.xlsx", sheet{sample});
    HD=data(:,27); HD1=data1(:,28); HD2=data2(:,29); HD3=data3(:,30);
    group1=[repmat({'twostate'}, size(HD,1),1);  repmat({'threestate'}, size(HD1,1),1);  repmat({'crosstalk'}, size(HD2,1),1);  repmat({'cross-talk three-state'}, size(HD3,1),1)];
    figure(2)
    boxplot([HD;  HD1;  HD2;  HD3], group1)
    title('HD')

    L = length(data(:,1));
     P0 = zeros(L, 8); Pm = nan(L, 8); m = nan(L, 8);
    for zushu = 1:L
        models = [1 2 3 4]; MM = [19 21 20 22]; CN = 10000;
      for ii = 1 : length(models)
        model = models(ii); MM1 = MM(ii);
        if model == 1 
           parameter = data(zushu,MM1:MM1+2);
           ydata = data01(zushu,:);
        elseif model == 2
            parameter = data2(zushu,MM1:MM1+2);
            ydata = data21(zushu,:);
        elseif model == 3 
            parameter = data1(zushu,MM1:MM1+2);
            ydata = data11(zushu,:);
        elseif model == 4
            parameter = data3(zushu,MM1:MM1+2);
            ydata = data31(zushu,:);
        end
        for i = 1:max(ydata)
           yyy(i) = sum(ydata==(i-1))/CN;
        end
        xdata=[]; tdata=[]; tt=[];
        xdata=yyy;
        tt=CN.*xdata;
        S=numel(xdata);
        for i=1:1:S;
            tdata(i)=i-1;
        end
        P0(zushu,ii) = xdata(1);
        dist = xdata;
        N = length(dist);
        for i=3:N-1
          if dist(i-1)<dist(i) && dist(i)>=dist(i+1)
             m(zushu,ii)=i-1;    
             Pm(zushu,ii)=dist(i);
           break
           end
        end
        clear dist yyy

        M = S-1;
        dist = fsp_model(parameter, 1, M);
        P0(zushu,ii+4) = dist(1);
        for i=3:N-1
          if dist(i-1)<dist(i) && dist(i)>=dist(i+1)
             m(zushu,ii+4)=i-1;    
             Pm(zushu,ii+4)=dist(i);
           break
           end
        end
        clear dist
      end
    end
    for i =1 :4
        delta_P0(i) = P0(:,i) - P0(:,i+4);
        delta_m(i) = m(:,i) - m(:,i+4);
        delta_Pm(i) = Pm(:,i) - Pm(:,i+4);
    end
    group2=[repmat({'twostate'}, size(clearnan(delta_P0(:,1)),1),1);  repmat({'threestate'}, size(clearnan(delta_P0(:,2)),1),1);  repmat({'crosstalk'}, size(clearnan(delta_P0(:,3)),1),1);  repmat({'cross-talk three-state'}, size(clearnan(delta_P0(:,4)),1),1)];
    figure(3)
    boxplot([clearnan(delta_P0(:,1));  clearnan(delta_P0(:,2));  clearnan(delta_P0(:,3));  clearnan(delta_P0(:,4))], group2)
    title('P0')
    clear group2
    group2=[repmat({'twostate'}, size(clearnan(delta_m(:,1)),1),1);  repmat({'threestate'}, size(clearnan(delta_m(:,2)),1),1);  repmat({'crosstalk'}, size(clearnan(delta_m(:,3)),1),1);  repmat({'cross-talk three-state'}, size(clearnan(delta_m(:,4)),1),1)];
    figure(4)
    boxplot([clearnan(delta_m(:,1));  clearnan(delta_m(:,2));  clearnan(delta_m(:,3));  clearnan(delta_m(:,4))], group2)
    title('m')
    clear group2
    group2=[repmat({'twostate'}, size(clearnan(delta_Pm(:,1)),1),1);  repmat({'threestate'}, size(clearnan(delta_Pm(:,2)),1),1);  repmat({'crosstalk'}, size(clearnan(delta_Pm(:,3)),1),1);  repmat({'cross-talk three-state'}, size(clearnan(delta_Pm(:,4)),1),1)];
    figure(5)
    boxplot([clearnan(delta_Pm(:,1));  clearnan(delta_Pm(:,2));  clearnan(delta_Pm(:,3));  clearnan(delta_Pm(:,4))], group2)
    title('P0')
    clear group2
end


%%%% figS3(b)
if sample <= 4
HD=[data(:,27) data(:,40) data(:,52)  data(:,66)];  %%two-state
HD1=[data1(:,28) data1(:,41) data1(:,53)  data1(:,67)];  %%three-state
HD2=[data2(:,29) data2(:,42) data2(:,54) data2(:,68)];  %%cross-talk
HD3=[data3(:,30) data3(:,43) data3(:,55)  data3(:,69)]; %%cross-talk three-state

AICc1=[data(:,22) data(:,35) data(:,47)    data(:,61)];
AICc2=[data1(:,23)   data1(:,36) data1(:,48)  data1(:,62)]; 
AICc3=[data2(:,24)  data2(:,37) data2(:,49)  data2(:,63)]; 
AICc4=[data3(:,25) data3(:,38) data3(:,50)    data3(:,64)]; 
AICc11=min(AICc1,[],2);
AICc21=min(AICc2,[],2);
AICc31=min(AICc3,[],2);
AICc41=min(AICc4,[],2);

HD1_tw=HD(AICc11== AICc1(:,1) | AICc11== AICc1(:,3), :); HD1_c=HD(AICc11== AICc1(:,2) | AICc11== AICc1(:,4), :);
HD2_tw=HD1(AICc21== AICc2(:,1) | AICc21== AICc2(:,3), :); HD2_c=HD1(AICc21== AICc2(:,2) | AICc21== AICc2(:,4), :);
HD3_tw=HD2(AICc31== AICc3(:,1) | AICc31== AICc3(:,3), :); HD3_c=HD2(AICc31== AICc3(:,2) | AICc31== AICc3(:,4), :);
HD4_tw=HD3(AICc41== AICc4(:,1) | AICc41== AICc4(:,3), :); HD4_c=HD3(AICc41== AICc4(:,2) | AICc41== AICc4(:,4), :);

r_HD1_tw=min([HD1_tw(:,1); HD1_tw(:,3)],[],2)./ min([HD1_tw(:,2); HD1_tw(:,4)],[],2);     r_HD1_c=min([HD1_c(:,1); HD1_c(:,3)],[],2)./ min([HD1_c(:,2); HD1_c(:,4)],[],2); 
r_HD2_tw=min([HD2_tw(:,1); HD2_tw(:,3)],[],2)./ min([HD2_tw(:,2); HD2_tw(:,4)],[],2);     r_HD2_c=min([HD2_c(:,1); HD2_c(:,3)],[],2)./ min([HD2_c(:,2); HD2_c(:,4)],[],2); 
r_HD3_tw=min([HD3_tw(:,1); HD3_tw(:,3)],[],2)./ min([HD3_tw(:,2); HD3_tw(:,4)],[],2);     r_HD3_c=min([HD3_c(:,1); HD3_c(:,3)],[],2)./ min([HD3_c(:,2); HD3_c(:,4)],[],2); 
r_HD4_tw=min([HD4_tw(:,1); HD4_tw(:,3)],[],2)./ min([HD4_tw(:,2); HD4_tw(:,4)],[],2);     r_HD4_c=min([HD4_c(:,1); HD4_c(:,3)],[],2)./ min([HD4_c(:,2); HD4_c(:,4)],[],2);  
group2 = [
    repmat({'twostate_tw'}, size(r_HD1_tw,1), 1);   
    repmat({'twostate_c'}, size(r_HD1_c,1), 1);    
    repmat({'threestate_tw'}, size(r_HD2_tw,1), 1); 
    repmat({'threestate_c'}, size(r_HD2_c,1), 1);   
    repmat({'crosstalk_tw'}, size(r_HD3_tw,1), 1); 
    repmat({'crosstalk_c'}, size(r_HD3_c,1), 1);   
    repmat({'crosstalk_threestate_tw'}, size(r_HD4_tw,1), 1); 
    repmat({'crosstalk_threestate_c'}, size(r_HD4_c,1), 1)    
];
positions = [0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2];
figure(2)
boxplot([r_HD1_tw; r_HD1_c; r_HD2_tw; r_HD2_c; r_HD3_tw; r_HD3_c; r_HD4_tw; r_HD4_c], group2, 'Positions', positions, 'Widths', 0.3);
hold on
kk=1.05;
plot([kk kk kk kk kk kk kk kk ])
title(sheet{sample})
ylim([0.97, 1.4])

end

%%%% figS3(a)
clc
clear all

model = 4; %% 1: two-state; 2: crosstalk; 3: threestate; 4:crosstalk-threestate
if model == 1
    filename = 'twostate_bimodal_fit.xlsx';
elseif model == 2
    filename = 'crosstalk_bimodal_fit.xlsx';
elseif model == 3
    filename = 'threestate_bimodal_fit.xlsx';
elseif model == 4
    filename = 'crosstalk_threestate_bimodal_fit.xlsx';
end

sheet = {'100', '250', '500', '1000'};
[data,name]=xlsread(filename, sheet{1});
[data1,name1]=xlsread(filename, sheet{2});
[data2,name2]=xlsread(filename, sheet{3});
[data3,name3]=xlsread(filename, sheet{4});

col = [22 35 47  61]; idx =[1 3 2 4]; i = idx(model) - 1;
AICc1=[data(:,col(1)+i)  data(:,col(2)+i)  data(:,col(3)+i)  data(:,col(4)+i)];
AICc2=[data1(:,col(1)+i)  data1(:,col(2)+i)  data1(:,col(3)+i)  data1(:,col(4)+i)]; 
AICc3=[data2(:,col(1)+i)  data2(:,col(2)+i)  data2(:,col(3)+i) data2(:,col(4)+i)]; 
AICc4=[data3(:,col(1)+i)  data3(:,col(2)+i)  data3(:,col(3)+i)  data3(:,col(4)+i)]; 
AICc11=min(AICc1,[],2);
AICc21=min(AICc2,[],2);
AICc31=min(AICc3,[],2);
AICc41=min(AICc4,[],2);
nn = length(data(:,1));
for i=1:4
    M100(i)=length(find(AICc11== AICc1(:,i)))/nn;
    M250(i)=length(find(AICc21== AICc2(:,i)))/nn;
    M500(i)=length(find(AICc31== AICc3(:,i)))/nn;
    M1000(i)=length(find(AICc41== AICc4(:,i)))/nn;
end

x=[1, 2,3 ,4];
y=[M1000; M500; M250; M100];

figure(4)
h=bar(x,y,1);
set(gca,'XTickLabel',{'1000', '500', '250', '100'})
xlabel('Sample size(N)')
ylim([0,1])
legend(' two-state model', 'cross-talk model', 'three-state model', 'cross-talk three-state model')
title('data')

