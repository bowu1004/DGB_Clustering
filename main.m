clear all;warning off;
%
% Copyright (c) 2017, Bo Wu, Bogdan M.Wilamowski
% All rights reserved.
%
% Project Title: A Fast Density and Grid Based Clustering Method for Data with Arbitrary Shapes and Noise
% Authors: Bo Wu, Bogdan M.Wilamowski
%
% Developer: Bo Wu, Bogdan M.Wilamowski
%
% Contact Info: bowu@auburn.edu, wilam@ieee.org
%
%

%% Choose input data and specify parameters
data_opt = 2;
switch data_opt
    case 1,
        inp_Data = load('./datasets/fig2_panelC.dat'); %inp_Data = load('fig2_panelB.dat');
        shouldbe_no_cluster = 5;
        no_grid=20;
        thre_grid_length=1.1; %search neighbor within distance <=1.1
        cutoff_factor=0.19;  %0.19*max_density  
        noise_thre = 2.5;
    case 2,
        inp_Data=load('./datasets/pathbased.txt');
        shouldbe_no_cluster = 3;
        no_grid=20;
        thre_grid_length=1.1;
        cutoff_factor=0.21;
        noise_thre = 2.5;
    case 3,
        inp_Data=load('./datasets/flame.txt');
        shouldbe_no_cluster = 2;
        no_grid=11; 
        thre_grid_length=1.1;
        cutoff_factor=0.42; 
        noise_thre = 0;
    case 4,
        inp_Data=load('./datasets/spiral.txt');
        shouldbe_no_cluster = 3;
        thre_grid_length=1.1;
        no_grid=22;
        cutoff_factor=0.07;
        noise_thre = 0;   
    case 5,
        inp_Data=load('./datasets/t4.8k.dat');
        result_title=sprintf('chameleon-data: t4.8k.dat');
        shouldbe_no_cluster = 6;
        no_grid=40;
        thre_grid_length=1.1;
        cutoff_factor=0.37; 
        noise_thre = 2.5;
    case 6,
        inp_Data=load('./datasets/t7.10k.dat');
        result_title=sprintf('chameleon-data: t7.10k.dat');
        shouldbe_no_cluster = 9;
        no_grid=40;
        thre_grid_length=1.1;
        cutoff_factor=0.4; 
        noise_thre = 2.5;
    case 7,
        inp_Data=load('./datasets/t8.8k.dat');
        result_title=sprintf('chameleon-data: t8.8k.dat');
        shouldbe_no_cluster = 8;
        no_grid=70;
        thre_grid_length=1.1; 
        cutoff_factor=0.16;
        noise_thre = 2.5;
end;
%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap is defined precisely, or use the default colormap
color_map={[1 0 0],[0 0.8 0.1],[0 0 1],[0.368 0.149 0.071],[1 0 1],[0 1 1],[0.8 0.05 0.35],...
    [0.43 0.83 0.68],[0.64 0.36 0.78],[0.6 0.2 0.7],[0.14 0.78 0.34],[0.3 0.4 0.5],...
    [0.5 0.7 0.3],[0.65 0.17 0.4],[1 1 0],...
    [148/255 148/255 148/255]};
no_Pattern = size(inp_Data,1);
x=inp_Data(:,1); xmin=min(x);xmax=max(x);
y=inp_Data(:,2); ymin=min(y);ymax=max(y);
figure(2),clf,plot(x,y,'.');title('Original Dataset');
%% Build the grid
%
% fuzzy clustering-I,different on the strategy of detection noise
%
stat_grid=[];
inp_all_mat=[];
dis_cnt=zeros(no_grid,no_grid);

dens_leftbottom=zeros(no_grid,no_grid);
dens_rightbottom=zeros(no_grid,no_grid);
dens_upperright=zeros(no_grid,no_grid);
dens_upperleft=zeros(no_grid,no_grid);

for i=1:no_Pattern
    x(i)=x(i)*(no_grid-1)/(xmax-xmin)+1-xmin*(no_grid-1)/(xmax-xmin); %start from 1 to no_Grid
    y(i)=y(i)*(no_grid-1)/(ymax-ymin)+1-ymin*(no_grid-1)/(ymax-ymin); %start from 1,because index>0
    
    x_minor=floor(x(i));
    x_larger=ceil(x(i));
    y_minor=floor(y(i));
    y_larger=ceil(y(i));
    
    if x_minor<1
        x_minor=1;
    end;
    if y_minor<1
        y_minor=1;
    end;
    if x_larger>no_grid
        x_larger=no_grid;
    end;
    if y_larger>no_grid
        y_larger=no_grid;
    end;
    
    dens_leftbottom(y_minor,x_minor)=dens_leftbottom(y_minor,x_minor)+(1-y(i)+y_minor)*(1-x(i)+x_minor);
    dens_rightbottom(y_minor,x_larger)=dens_rightbottom(y_minor,x_larger)+(1-y(i)+y_minor)*(1-x_larger+x(i));
    dens_upperright(y_larger,x_larger)=dens_upperright(y_larger,x_larger)+(1-y_larger+y(i))*(1-x_larger+x(i));
    dens_upperleft(y_larger,x_minor)=dens_upperleft(y_larger,x_minor)+(1-y_larger+y(i))*(1-x(i)+x_minor);
end;

dis_cnt=dens_leftbottom+dens_rightbottom+dens_upperright+dens_upperleft; %1+2+3+4 corners

figure(4),clf,plot(x,y,'.');grid on;hold on;
X=0:1:no_grid;
Y=0:1:no_grid;
set(gca,'xtick',X,'ytick',Y);
s=sprintf('Mapping into %d*%d grid',no_grid,no_grid);
title(s);
%
% Show density map
%
range=[1 no_grid 1 no_grid];
step_x=1;
step_y=1;
[xx,yy]=meshgrid(range(1):step_x:range(2),range(3):step_y:range(4));
figure(44),clf,mesh(xx,yy,dis_cnt);title('density map');
%% Statistical resutls analysis
dis_cnt_ori=dis_cnt;
% dis_cnt_ori=dis_cnt_update;

NN=no_grid*no_grid;
cl=zeros(1,NN); %cluster index
cl_mat=zeros(no_grid,no_grid);
no_cluster=1;
marked_number=[1];

rho=reshape(dis_cnt,1,no_grid*no_grid); %store column followed by column
[rho_val, rho_ord]=sort(rho,'descend'); %sort density rho
MAX_DENSITY=max(rho); %maxium counts
maxd=MAX_DENSITY;

    
edge_factor_mat=[];
tic;
% find the mountain ridges, or the number of clusters automatically
while  (maxd>MAX_DENSITY*cutoff_factor)    
    cluster_add_flag=1; % 1:add a new cluster;0:don't add new cluster.
    len_mk_no=length(marked_number);
    edge_factor=cutoff_factor*maxd;
    edge_factor_mat=[edge_factor_mat maxd];
    for i=1:len_mk_no     
        % %-------------------get_matrix_position-----------------% %
        cur_col_i=floor(rho_ord(i)/no_grid)+1; %x-coordinate
        cur_row_i=mod(rho_ord(i),no_grid); %y-coordinate
        if cur_row_i==0
            cur_col_i=cur_col_i-1;
            cur_row_i=no_grid;
        end;
        % %-------------------------------------------------------% %
        
        cl_mat(cur_row_i,cur_col_i)=no_cluster;
        
        for j=(length(marked_number)+1):NN
            if rho_val(j)>edge_factor %thre_rho,edge detection
                % %-------------------get_matrix_position-------------% %
                cur_col_j=floor(rho_ord(j)/no_grid)+1; %x-coordinate
                cur_row_j=mod(rho_ord(j),no_grid); %y-coordinate
                if cur_row_j==0
                    cur_col_j=cur_col_j-1;
                    cur_row_j=no_grid;
                end;
                % %---------------------------------------------------% %
                if cl_mat(cur_row_j,cur_col_j)==0 %haven't marked yet
                    % check if it's neighbor or not
                    dis_tmp=abs(cur_row_i-cur_row_j)+abs(cur_col_i-cur_col_j);
                    if dis_tmp<thre_grid_length
                        cl_mat(cur_row_j,cur_col_j)=cl_mat(cur_row_i,cur_col_i);
                        dis_cnt(cur_row_j,cur_col_j)=maxd;
                        marked_number=[marked_number j]; %marked grids
                        cluster_add_flag=0;
                    end;
                end;
            end;
        end;
    end;
    
    if (cluster_add_flag>0)
        % add a new cluster
        for ii=1:length(marked_number)
            % %-------------------get_matrix_position-----------------% %
            cur_col_ii=floor(rho_ord(ii)/no_grid)+1; %x-coordinate
            cur_row_ii=mod(rho_ord(ii),no_grid); %y-coordinate
            if cur_row_ii==0
                cur_col_ii=cur_col_ii-1;
                cur_row_ii=no_grid;
            end;
            % %-------------------------------------------------------% %
            dis_cnt(cur_row_ii,cur_col_ii)=0;
        end;
        no_cluster=no_cluster+1;
        marked_number=[1];
    end;
    
    rho=reshape(dis_cnt,1,no_grid*no_grid); %store column followed by column
    [rho_val, rho_ord]=sort(rho,'descend'); %sort density rho
    maxd=max(rho); %maxium counts
end;
% show results
figure(22),clf,mesh(xx,yy,cl_mat);title('Mountain Ridges');colorbar;

no_cluster=no_cluster-1;
disp(['number of clusters is: ' num2str(no_cluster)]);
%% Noise filtering
figure(34),clf;
for i=1:no_Pattern
    node_i_row=round(y(i));
    node_i_col=round(x(i));
    tmp_cl=cl_mat(node_i_row,node_i_col);
    if tmp_cl==0
        tmp_cl=16;
        if  dis_cnt_ori(node_i_row,node_i_col)>noise_thre
            server_node_val_0=0;
            non_zero_nei_node=1;
            x_minor=floor(x(i));
            x_larger=ceil(x(i));
            y_minor=floor(y(i));
            y_larger=ceil(y(i));
            
            if x_minor<1
                x_minor=1;
            end;
            if y_minor<1
                y_minor=1;
            end;
            if x_larger>no_grid
                x_larger=no_grid;
            end;
            if y_larger>no_grid
                y_larger=no_grid;
            end;
            % four corner labels
            % %-----------------4     3-----------------% %
            % %--------------sparse pattern-------------% %
            % %-----------------1     2-----------------% %
            if cl_mat(y_minor,x_minor)>0
                tmp_cl=cl_mat(y_minor,x_minor);
                corner_position=1;
                server_node_val_0=server_node_val_0+dens_leftbottom(y_minor,x_minor);
                non_zero_nei_node=non_zero_nei_node+1;
            end;
            if cl_mat(y_minor,x_larger)>0
                tmp_cl=cl_mat(y_minor,x_larger);
                corner_position=2;
                server_node_val_0=server_node_val_0+dens_rightbottom(y_minor,x_larger);
                non_zero_nei_node=non_zero_nei_node+1;
            end;
            if cl_mat(y_larger,x_larger)>0
                tmp_cl=cl_mat(y_larger,x_larger);
                corner_position=3;
                server_node_val_0=server_node_val_0+dens_upperright(y_larger,x_larger);
                non_zero_nei_node=non_zero_nei_node+1;
            end;
            if cl_mat(y_larger,x_minor)>0
                tmp_cl=cl_mat(y_larger,x_minor);
                corner_position=4;
                server_node_val_0=server_node_val_0+dens_upperleft(y_larger,x_minor);
                non_zero_nei_node=non_zero_nei_node+1;
            end;
            % !!!
            if non_zero_nei_node==2 %only one of four neighbor nodes has a label
                if (node_i_row==y_minor) && (node_i_col==x_minor) %the nearest node, corner 1
                    client_node_val=dens_leftbottom(y_minor,x_minor);
                elseif (node_i_row==y_minor) && (node_i_col==x_larger) %the nearest node, corner 2
                    client_node_val=dens_rightbottom(y_minor,x_larger);
                elseif (node_i_row==y_larger) && (node_i_col==x_larger) %the nearest node, corner 3
                    client_node_val=dens_upperright(y_larger,x_larger);
                else %the nearest node, corner 4
                    client_node_val=dens_upperleft(y_larger,x_minor);
                end;
                if client_node_val<noise_thre
                    tmp_cl=16; %noise
                end;
            end;
            if non_zero_nei_node>2
                if server_node_val_0<noise_thre
                    tmp_cl=16; %noise
                end;
            end;
        end;
    end;
    plot(inp_Data(i,1),inp_Data(i,2),'.','color',cell2mat(color_map(tmp_cl)),'markersize',10);hold on;
end;
disp(['processing time: ' num2str(toc) 'sec']);