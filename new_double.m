%双摆
clear
clc
close all
%输入
N=2;%双摆
m1=1;
m2=1;
l1=1;
l2=1;
g=9.8;
Input=[N,m1,m2,l1,l2,g];
%初始条件和时间设置
y0=[pi/2;pi/2;0;0];%这里全部是弧度值。分别代表[摆1与垂面夹角，摆2与垂面夹角，摆1角动量，摆2角动量]
% h=1e-2;
h=0.0001;
x0=0:h:20;
%代入到ODE求解器中
[y1,Output]=ODE_RK4_hyh(x0,h,y0,Input);

%提取出角度
tN=size(y1,2);
th1=y1(1,:);
th2=y1(2,:);
p1=y1(3,:);
p2=y1(4,:);
% 计算导数值
dth1 = diff(th1);
dth2 = diff(th2);
dp1 = diff(p1);
dp2 = diff(p2);

%计算出关节坐标
CX1_A=zeros(1,tN);
CX1_B=CX1_A+l1*sin(th1);
CY1_A=zeros(1,tN);
CY1_B=CY1_A-l1*cos(th1);

CX2_A=CX1_B;
CX2_B=CX2_A+l2*sin(th2);
CY2_A=CY1_B;
CY2_B=CY2_A-l2*cos(th2);


%绘图
n=1;
figure()
set(gcf,'position',[488   342   400   300])
for k=1:4:length(x0) %这里4步一显示时间帧
    clf
    xlim([-2,2])
    ylim([-2,2])
    hold on
    %绘制摆
    plot([CX1_A(k),CX1_B(k)],[CY1_A(k),CY1_B(k)],'color','k','LineWidth',1.5)
    plot([CX2_A(k),CX2_B(k)],[CY2_A(k),CY2_B(k)],'color','k','LineWidth',1.5)
    %绘制轨线
    if k>200
        n=n+1;
    end
    Nm=k-n+1;
    %轨迹1
    F_color=[1,0,0];
    F_color=F_color*0.6+[1,1,1]*0.4*0.999;
    cdata=[linspace(1,F_color(1),Nm+1)',linspace(1,F_color(2),Nm+1)',linspace(1,F_color(3),Nm+1)'];
    cdata=reshape(cdata,Nm+1,1,3);
    if k>3
        patch([CX1_B(n:k),NaN],[CY1_B(n:k),NaN],1:Nm+1,'EdgeColor','interp','Marker','none',...
          'MarkerFaceColor','flat','CData',cdata,'LineWidth',1.5);
    end
    %轨迹2
    F_color=[0,0,1];
    F_color=F_color*0.6+[1,1,1]*0.4*0.999;
    cdata=[linspace(1,F_color(1),Nm+1)',linspace(1,F_color(2),Nm+1)',linspace(1,F_color(3),Nm+1)'];
    cdata=reshape(cdata,Nm+1,1,3);
    if k>3
        patch([CX2_B(n:k),NaN],[CY2_B(n:k),NaN],1:Nm+1,'EdgeColor','interp','Marker','none',...
          'MarkerFaceColor','flat','CData',cdata,'LineWidth',1.5);
    end
    hold off
    pause(0.05)
    %可以在这里添加输出动图的程序
end

function [y,Output]=ODE_RK4_hyh(x,h,y0,Input)
%4阶RK方法
%h间隔为常数的算法
y=zeros(size(y0,1),size(x,2));
y(:,1)=y0;
for ii=1:length(x)-1
    yn=y(:,ii);
    xn=x(ii);
    K1=Fdydx(xn    ,yn       ,Input);
    K2=Fdydx(xn+h/2,yn+h/2*K1,Input);
    K3=Fdydx(xn+h/2,yn+h/2*K2,Input);
    K4=Fdydx(xn+h  ,yn+h*K3  ,Input);
    y(:,ii+1)=yn+h/6*(K1+2*K2+2*K3+K4);
end
Output=[];
end

function dydx=Fdydx(x,y,Input)
%将原方程整理为dy/dx=F(y,x)的形式
%输入Input整理
m1=Input(2);
m2=Input(3);
l1=Input(4);
l2=Input(5);
g=Input(6);
%输入
theta_1=y(1);%角度1
theta_2=y(2);%角度2
p_1=y(3);%角动量1
p_2=y(4);%角动量2
%利用拉格朗日法得到的方程
dtheta_1 = (12*(l2*p_1 + 3*l2*p_1*cos(theta_2)^2 + 3*l2*p_1*sin(theta_2)^2 - 6*l1*p_2*sin(theta_1)*sin(theta_2) - 6*l1*p_2*cos(theta_1)*cos(theta_2)))/(l1^2*(l2*m1 + 3*l2*m1*cos(theta_1)^2 + 3*l2*m1*cos(theta_2)^2 + 12*l2*m2*cos(theta_1)^2 + 3*l2*m1*sin(theta_1)^2 + 3*l2*m1*sin(theta_2)^2 + 12*l2*m2*sin(theta_1)^2 + 9*l2*m1*cos(theta_1)^2*cos(theta_2)^2 + 9*l2*m1*cos(theta_1)^2*sin(theta_2)^2 + 9*l2*m1*cos(theta_2)^2*sin(theta_1)^2 + 36*l2*m2*cos(theta_1)^2*sin(theta_2)^2 + 36*l2*m2*cos(theta_2)^2*sin(theta_1)^2 + 9*l2*m1*sin(theta_1)^2*sin(theta_2)^2 - 72*l2*m2*cos(theta_1)*cos(theta_2)*sin(theta_1)*sin(theta_2)));

dtheta_2 =(12*(l1*m1*p_2 + 3*l1*m1*p_2*cos(theta_1)^2 + 12*l1*m2*p_2*cos(theta_1)^2 + 3*l1*m1*p_2*sin(theta_1)^2 + 12*l1*m2*p_2*sin(theta_1)^2 - 6*l2*m2*p_1*cos(theta_1)*cos(theta_2) - 6*l2*m2*p_1*sin(theta_1)*sin(theta_2)))/(l1*m2*(l2^2*m1 + 3*l2^2*m1*cos(theta_1)^2 + 3*l2^2*m1*cos(theta_2)^2 + 12*l2^2*m2*cos(theta_1)^2 + 3*l2^2*m1*sin(theta_1)^2 + 3*l2^2*m1*sin(theta_2)^2 + 12*l2^2*m2*sin(theta_1)^2 + 9*l2^2*m1*cos(theta_1)^2*cos(theta_2)^2 + 9*l2^2*m1*cos(theta_1)^2*sin(theta_2)^2 + 9*l2^2*m1*cos(theta_2)^2*sin(theta_1)^2 + 36*l2^2*m2*cos(theta_1)^2*sin(theta_2)^2 + 36*l2^2*m2*cos(theta_2)^2*sin(theta_1)^2 + 9*l2^2*m1*sin(theta_1)^2*sin(theta_2)^2 - 72*l2^2*m2*cos(theta_1)*cos(theta_2)*sin(theta_1)*sin(theta_2)));

dp1=- (m2*(2*dtheta_1*l1*sin(theta_1)*(dtheta_1*l1*cos(theta_1) + (dtheta_2*l2*cos(theta_2))/2) - 2*dtheta_1*l1*cos(theta_1)*(dtheta_1*l1*sin(theta_1) + (dtheta_2*l2*sin(theta_2))/2)))/2 - (3*g*l1*m1*sin(theta_1))/2;
dp2=- (m2*(dtheta_2*l2*sin(theta_2)*(dtheta_1*l1*cos(theta_1) + (dtheta_2*l2*cos(theta_2))/2) - dtheta_2*l2*cos(theta_2)*(dtheta_1*l1*sin(theta_1) + (dtheta_2*l2*sin(theta_2))/2)))/2 - (g*l2*m1*sin(theta_2))/2;
%整理输出
dydx=zeros(4,1);
dydx(1)=dtheta_1;
dydx(2)=dtheta_2;
dydx(3)=dp1;
dydx(4)=dp2;
end

