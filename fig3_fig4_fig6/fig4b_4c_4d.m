clc
clear all

kk = 1; % 1: bimodal; 2: unimodal

[data1,name1]=xlsread("ESM_kd.xls");
kd=data1(:,7);      
Gname1=name1(2:end,4);

[data,name]=xlsread("41586_2024_7517_MOESM4_ESM.xlsx");
index=find(startsWith(name(:,1), 'GN-'));
bf=data(index-1,7);                       %%% nascent bf
bz=data(index-1,8);
Gname=extractAfter(name(index,1),'GN-');  %%% nascent gene name
clear data name data1 name1

[common_genes, idx1, idx2] = intersect(Gname1, Gname, 'stable');
gene_kd = kd(idx1);
gene_bf = bf(idx2);
gene_bz = bz(idx2);
clear data1 name1 data name index kd bf bz Gname1 Gname idx1 idx2

[data,name]=xlsread("Embryonic_c57_fitmethod_nofixp.csv");
[~,name1] = xlsread("Supplement Table S1.xlsx", 'E c57');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
name21 = name(2:end,1);
if kk == 1
    data2=data(row,:);  name2=name21(row,1);
elseif kk == 2
    data(row,:) = []; name21 = name(2:end,1); name21(row,:) =[];  
    data2 = data;  name2 = name21;
end
clear data name   name21

[~, idx1, idx2] = intersect(common_genes, name2, 'stable');
data_E1=data2(idx2,:);
kd_E1 = gene_kd(idx1);
bf_E1=gene_bf(idx1); 
bz_E1=gene_bz(idx1);
clear  idx1 idx2 data2 name2  

[data,name]=xlsread("Embryonic_cast_fitmethod_nofixp.csv");
[~,name1] = xlsread("Supplement Table S1.xlsx", 'E cast');
[~, row, ~] =intersect(name(2:end,1), name1(2:end,1), 'stable');
name21 = name(2:end,1);
if kk == 1
    data2=data(row,:);  name2=name21(row,1);
elseif kk == 2
    data(row,:) = [];  name21(row,:) =[];  
    data2 = data;  name2 = name21;
end
clear data name  idx1 idx2 kd name21

[~, idx1, idx2] = intersect(common_genes, name2, 'stable');
data_E2=data2(idx2,:);
kd_E2 = gene_kd(idx1);
bf_E2=gene_bf(idx1); 
bz_E2=gene_bz(idx1);
clear common_genes idx1 idx2 data2 name2 gene_bz gene_bf gene_kd row

AICc_E1=[data_E1(:,10) data_E1(:,24)];    AICc_E11=min(AICc_E1,[],2);
row1= find(AICc_E1(:,1)==AICc_E11 ); 
row11 = find(AICc_E1(:,2)==AICc_E11 ); 
para_E1 = [ data_E1(:, 6:8), data_E1(:, 18:22) ];

AICc_E2=[data_E2(:,10) data_E2(:,24)];    AICc_E21=min(AICc_E2,[],2);
row2= find(AICc_E2(:,1)==AICc_E21); 
row21 = find(AICc_E2(:,2)==AICc_E21); 
para_E2 = [ data_E2(:, 6:8), data_E2( :, 18:22) ];

bf_e_E1 = kd_E1 .* [ para_E1(:, 1 ),  1./( para_E1(:,6)./para_E1(:,4) + (1-para_E1(:,6))./para_E1(:,5) )  ];
bz_e_E1 = [para_E1(: ,3)./para_E1(: ,2),  para_E1(: ,8)./para_E1(:,7) ];
bf_e_E2 = kd_E2 .* [ para_E2(: , 1 ),  1./( para_E2(: ,6)./para_E2(: ,4) + (1-para_E2(:,6))./para_E2(:,5) )  ];
bz_e_E2 = [para_E2(:,3)./para_E2(: ,2),  para_E2(: ,8)./para_E2(: ,7) ];

r_bf_E1 = abs(bf_e_E1 - bf_E1);
r_bf_E2 = abs(bf_e_E2 - bf_E2);
r_bz_E1 = abs(bz_e_E1 - bz_E1);
r_bz_E2 = abs(bz_e_E2 - bz_E2);

if kk == 1
%%% fig4b
figure(1)
group1 = [
    repmat({'c57'}, size(r_bf_E1(row11,1),1), 1);   
    repmat({'c57c'}, size(r_bf_E1(row11,2), 1), 1); 
    repmat({'cast'}, size(r_bf_E2(row21,1),1), 1);  
    repmat({'castc'}, size(r_bf_E2(row21,2),1), 1)  
];
positions = [0.8, 1.2, 1.8, 2.2]; 
h = boxplot([r_bf_E1(row11,1); r_bf_E1(row11,2); r_bf_E2(row21,1); r_bf_E2(row21,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0, 5])
ylable('bf')
xlabel('mean')

group2 = [
    repmat({'c57'}, size(r_bz_E1(row11,1),1), 1);   
    repmat({'c57c'}, size(r_bz_E1(row11,2), 1), 1); 
    repmat({'cast'}, size(r_bz_E2(row21,1),1), 1);  
    repmat({'castc'}, size(r_bz_E2(row21,2),1), 1)  
];
figure(2)
positions = [0.8, 1.2, 1.8, 2.2]; 
h = boxplot([r_bz_E1(row11,1); r_bz_E1(row11,2); r_bz_E2(row21,1); r_bz_E2(row21,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0, 115])

clear group1 group2
%%% fig4c
figure(3)
group1 = [
    repmat({'c57'}, size(r_bf_E1(row1,1),1), 1);   
    repmat({'c57c'}, size(r_bf_E1(row1,2), 1), 1); 
    repmat({'cast'}, size(r_bf_E2(row2,1),1), 1);  
    repmat({'castc'}, size(r_bf_E2(row2,2),1), 1)  
];
positions = [0.8, 1.2, 1.8, 2.2]; 
h = boxplot([r_bf_E1(row1,1); r_bf_E1(row1,2); r_bf_E2(row2,1); r_bf_E2(row2,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0.5, 5])
ylable('bf')
xlabel('mean')

group2 = [
    repmat({'c57'}, size(r_bz_E1(row1,1),1), 1);   
    repmat({'c57c'}, size(r_bz_E1(row1,2), 1), 1); 
    repmat({'cast'}, size(r_bz_E2(row2,1),1), 1);  
    repmat({'castc'}, size(r_bz_E2(row2,2),1), 1)  
];
figure(4)
positions = [0.8, 1.2, 1.8, 2.2]; 
h = boxplot([r_bz_E1(row1,1); r_bz_E1(row1,2); r_bz_E2(row2,1); r_bz_E2(row2,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0, 185])
end

%%% fig4d
if kk == 2
group1 = [
    repmat({'c57'}, size(r_bf_E1(:,1),1), 1);   
    repmat({'c57c'}, size(r_bf_E1(:,2), 1), 1); 
    repmat({'cast'}, size(r_bf_E2(:,1),1), 1);  
    repmat({'castc'}, size(r_bf_E2(:,2),1), 1)  
];
figure(1)
positions = [0.8, 1.2, 1.8, 2.2]; 
h = boxplot([r_bf_E1(:,1); r_bf_E1(:,2); r_bf_E2(:,1); r_bf_E2(:,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0, 4.5])

group2 = [
    repmat({'c57'}, size(r_bz_E1(:,1),1), 1);    
    repmat({'c57c'}, size(r_bz_E1(:,2), 1), 1);  
    repmat({'cast'}, size(r_bz_E2(:,1),1), 1);   
    repmat({'castc'}, size(r_bz_E2(:,2),1), 1)  
];
figure(2)
positions = [0.8, 1.2, 1.8, 2.2];  
h1 = boxplot([r_bz_E1(:,1); r_bz_E1(:,2); r_bz_E2(:,1); r_bz_E2(:,2)], group1, 'Positions', positions,'Widths', 0.3);
ylim([0, 50])

end


