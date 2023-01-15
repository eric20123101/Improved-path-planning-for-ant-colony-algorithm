%% Improved path planning for ant colony algorithm

% author : Eric_lee

% e-mail : 282228095@qq.com

% date   : 2019/4/22

%% 构建栅格地图
clear;clc
G=rand(30)>0.8; % G 地形图为01矩阵，如果为1，表示障碍物
save('D:\MATLAB\R2016b\work\AGV\map.mat','G');%将生成的随机地图以mat格式保存到本地，便于和标准算法在同一个地图下进行对比，这个路径可以改为自己电脑的路径
%load map30.mat
h=rot90(abs(peaks(30)));
global MM
global Dir 
global Lgrid 
%Lgrid = input('请输入栅格粒径：');
Lgrid = 1;
MM=size(G,1); %MM为矩阵维数
figure(1);   
for i=1:MM
  for j=1:MM
  x1=(j-1)*Lgrid;y1=(MM-i)*Lgrid; 
  
  x2=j*Lgrid;y2=(MM-i)*Lgrid; 
  x3=j*Lgrid;y3=(MM-i+1)*Lgrid; 
  x4=(j-1)*Lgrid;y4=(MM-i+1)*Lgrid; 
  f=(max(max(h))-h(i,j))/max(max(h));
    if G(i,j)==1 
        fill([x1,x2,x3,x4],[y1,y2,y3,y4],[0.2,0.2,0.2]); hold on %栅格为1，填充为黑色
    else 
        fill([x1,x2,x3,x4],[y1,y2,y3,y4],[f,1,f]); hold on %栅格为0，填充为白色
    end 
  end 
end
axis([0,MM*Lgrid,0,MM*Lgrid]) 
grid on
%% 初始化地图信息
%{
Xinitial = input('请输入初始点的X坐标：');
Yinitial = input('请输入初始点的Y坐标：');
%}
Xinitial = 0.6;
Yinitial = 0.2;
[initial,ij_initial]= modify(Xinitial,Yinitial);
if max(ij_initial)>MM||G(ij_initial(1),ij_initial(2))==1
    error('初始点不能设在障碍物上或超出范围');
end
%{
Xdestination = input('请输入目标点的X坐标：');
Ydestination = input('请输入目标点的Y坐标：');
    %}
Xdestination = 29.4;
Ydestination = 29.3;
[destination,ij_destination]= modify(Xdestination,Ydestination);
if max(ij_destination)>MM||G(ij_destination(1),ij_destination(2))==1
    error('目标点不能设在障碍物上或超出范围');
end
%% 计算距离启发矩阵dis
dis = zeros(MM,MM);
for i=1:MM
  for j=1:MM
   x = (j-0.5)*Lgrid;
   y = (MM-i+0.5)*Lgrid;
   dis(i,j) = sqrt(sum(([x y]-destination).^2));
  end
end
%% 计算距离转移矩阵D
D=zeros(MM^2,8);   %行号表示栅格标号，列号表示邻接的8个方向的栅格号
Dir = [-MM-1,-1,MM-1,MM,MM+1,1,1-MM,-MM];
 for i = 1:MM^2     %8方向转移距离矩阵初步构建
     Dirn = Dir+i;
     if G(i)==1
             D(i,:)=inf;
             continue
     end
         for j = 1:8
             if  Dirn(j)<=0||Dirn(j)>MM^2        %出界的情况，暂且为0
                 continue 
             end
             if G(Dirn(j))==1
                 D(i,j) = inf;
             elseif mod(j,2)==0         %偶数方向为上下左右方向
                 D(i,j) = 1;
             elseif j==1 %左上方向的情况，保证路线不会擦障碍物边沿走过
                 if (G(Dirn(2))+G(Dirn(8))==0)
                   D(i,j) = 1.4; 
                 else
                   D(i,j) = inf;   
                 end
             elseif (Dirn(j-1)<=0||Dirn(j-1)>MM^2)||(Dirn(j+1)<=0||Dirn(j+1)>MM^2)%排除掉垂直方向的栅格出界的情况
                 continue
             elseif G(Dirn(j-1))+G(Dirn(j+1))==0    %其余三个斜方向
                 D(i,j) = 1.4;
             else
                 D(i,j) = inf;
             end
         end
     
 end
%% 创造边界
 num = 1:MM^2;
 obs_up = find(mod(num,MM)==1);
 obs_up = obs_up(2:end-1);
 D(obs_up,[1,2,3])=inf;
 obs_down = find(mod(num,MM)==0);
 obs_down = obs_down(2:end-1);
 D(obs_down,[5,6,7])=inf;
 D(2:MM-1,[1,7,8]) = inf;
 D(MM^2-MM+2:MM^2-1,[3,4,5])=inf;
 D(1,[1,2,3,7,8])=inf;
 D(MM,[1,5,6,7,8])=inf;
 D(MM^2-MM+1,[1,2,3,4,5])=inf;
 D(MM^2,[3,4,5,6,7])=inf;
%% 参数初始化
tic
NC_max=30; m=50;  Rho=0.3; Q=100; Omega=10; Mu=1;  u=10; Tau_min=10; Tau_max=40; Rho_min=0.2;
%% 绘制找到的最优路径
[R_best,F_best,L_best,T_best,S_best,S_ave,Shortest_Route,Shortest_Length]=improved(D,initial,destination,dis,h,NC_max,m,Rho,Omega,Mu,Q,u,Tau_min,Tau_max,Rho_min); %函数调用
%绘制找到的最优路径
j = ceil(Shortest_Route/MM);
i = mod(Shortest_Route,MM);
i(i==0) = MM;
x = (j-0.5)*Lgrid;
y = (MM-i+0.5)*Lgrid;
x = [initial(1) x destination(1)];
y = [initial(2) y destination(2)];
figure(1);
plot(x,y,'-r');
xlabel('x'); ylabel('y'); title('最佳路径');
grid on 
hold on
toc  %计算运行时间
%% 绘制收敛曲线
figure(2); iter=1:length(L_best);
plot(iter,L_best,'-r','LineWidth',1)
xlabel('迭代次数'); ylabel('各代最佳路线的长度');
axis([0,NC_max,25,90]);
grid on;hold on
figure(3); iter=1:length(L_best);
plot(iter,F_best*100,'-r','LineWidth',1)
xlabel('迭代次数'); ylabel('各代最佳路线的高度均方差*100'); 
axis([0,NC_max,0,30]);
grid on;hold on
figure(4); iter=1:length(L_best);
plot(iter,T_best,'-r','LineWidth',1)
xlabel('迭代次数'); ylabel('各代最佳路线的转弯次数'); 
axis([0,NC_max,5,50]);
grid on;hold on
figure(5); iter=1:length(L_best);
plot(iter,S_best,'r',iter,S_ave,'b');
xlabel('迭代次数');ylabel('各代最佳路线的综合指标及平均综合指标');
title('收敛性分析曲线')
axis([0,NC_max,70,200]);
grid on;hold on
toc


