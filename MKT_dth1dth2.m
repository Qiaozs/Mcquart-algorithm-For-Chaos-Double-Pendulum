% 将两个动力学方程同时引入构造目标函数，可以成功辨识，目标函数值收敛到0，双摆系统存在多解所以可能不能得到真实参数。
clc;
% 数据量
num = 2000;
% 单位阵   
I = eye(4);
% 迭代次数
N = 10000;
% 初始化
A = zeros(4,4);
Q = zeros(4,1);
a = zeros(4,4);
% 实际数据
x = [th1;th2;p1;p2];
y = [diff(th1)*10000;diff(th2)*10000];
%% 

% 初始化目标函数
J_aim = zeros(1,N);
v = zeros(4,N);

% 给定第一个参数集
parameter = zeros(4,N);
parameter(1,1) = 5;
parameter(2,1) = 1;
parameter(3,1) = 1;
parameter(4,1) = 1;
% 终止条件
e = 0.00000001;
    % 初始化计算数据
    dtheta_1 = zeros(1,num);
    dtheta_2 = zeros(1,num);
    dp_1 = zeros(1,num);
    dp_2 = zeros(1,num);
    % 初始化detla值
    detlay = zeros(2,num);
    compute_output = zeros(2,num);

% 计算第一个目标函数J_aim0
for j = 1 : num
    % 参数和数据传入
   m1 = parameter(1,1); m2 = parameter(2,1); l1 = parameter(3,1); l2 = parameter(4,1);
   theta_1 = x(1,j); theta_2 = x(2,j); p_1 = x(3,j); p_2 = x(4,j);
   dtheta_1(j)= (12*(2*l2*p_1 - 3*l1*p_2*cos(theta_1 - theta_2)))/(l1^2*l2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
   dtheta_2(j) = (12*(2*l1*m1*p_2 + 6*l1*m2*p_2 - 3*l2*m2*p_1*cos(theta_1 - theta_2)))/(l1*l2^2*m2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
   compute_output(:,j) = [dtheta_1(j);dtheta_2(j)];

   % 计算在当前参数下实际输出和辨识输出差
   detlay(:,j) = y(:,j)-compute_output(:,j); 
   % 根据当前参数计算目标函数值
   J_aim(1) = J_aim(1) +sum(power(detlay(:,j),2));
end

%% 

lambda = 1;
t = 2;
while true
    % 计算A和Q
    m1 = parameter(1,t-1); m2 = parameter(2,t-1); l1 = parameter(3,t-1); l2 = parameter(4,t-1);
    A = zeros(4,4);
    Q = zeros(4,1);
    for i = 1:num
        %数据传入
      theta_1 = x(1,i); theta_2 = x(2,i); p_1 = x(3,i); p_2 = x(4,i);

      dtheta_1(i) = (12*(2*l2*p_1 - 3*l1*p_2*cos(theta_1 - theta_2)))/(l1^2*l2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
     dtheta_2(i) = (12*(2*l1*m1*p_2 + 6*l1*m2*p_2 - 3*l2*m2*p_1*cos(theta_1 - theta_2)))/(l1*l2^2*m2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
   compute_output(:,i) = [dtheta_1(i);dtheta_2(i)];

   % 雅可比矩阵,用于计算预测函数的一阶微分，大小4×1
ddtheta1dm1 =-(96*(4*l2*p_1 - 6*l1*p_2*cos(theta_1 - theta_2)))/(l1^2*l2*(480*m1*m2 - 540*m2^2*cos(2*theta_1 - 2*theta_2) + 81*m2^2*cos(4*theta_1 - 4*theta_2) + 128*m1^2 + 531*m2^2 - 288*m1*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta1dm2 =-(48*(15*l2*p_1 - (63*l1*p_2*cos(theta_1 - theta_2))/4 - 9*l2*p_1*cos(2*theta_1 - 2*theta_2) + (27*l1*p_2*cos(3*theta_1 - 3*theta_2))/4))/(l1^2*l2*(480*m1*m2 - 540*m2^2*cos(2*theta_1 - 2*theta_2) + 81*m2^2*cos(4*theta_1 - 4*theta_2) + 128*m1^2 + 531*m2^2 - 288*m1*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta1dl1 =-(48*l2*p_1 - 36*l1*p_2*cos(theta_1 - theta_2))/(l1^3*l2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta1dl2 =(36*p_2*cos(theta_1 - theta_2))/(l1*l2^2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));  


ddtheta2dm1 =-(432*l1*p_2 - 576*l2*p_1*cos(theta_1 - theta_2) + 432*l1*p_2*cos(2*theta_1 - 2*theta_2))/(l1*l2^2*(480*m1*m2 - 540*m2^2*cos(2*theta_1 - 2*theta_2) + 81*m2^2*cos(4*theta_1 - 4*theta_2) + 128*m1^2 + 531*m2^2 - 288*m1*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta2dm2 =-(384*l1*m1^2*p_2 + 2160*l1*m2^2*p_2 - 756*l2*m2^2*p_1*cos(theta_1 - theta_2) + 1440*l1*m1*m2*p_2 - 1296*l1*m2^2*p_2*cos(2*theta_1 - 2*theta_2) + 324*l2*m2^2*p_1*cos(3*theta_1 - 3*theta_2) - 864*l1*m1*m2*p_2*cos(2*theta_1 - 2*theta_2))/(l1*l2^2*m2^2*(480*m1*m2 - 540*m2^2*cos(2*theta_1 - 2*theta_2) + 81*m2^2*cos(4*theta_1 - 4*theta_2) + 128*m1^2 + 531*m2^2 - 288*m1*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta2dl1 =(36*p_1*cos(theta_1 - theta_2))/(l1^2*l2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
ddtheta2dl2 =-(48*l1*m1*p_2 + 144*l1*m2*p_2 - 36*l2*m2*p_1*cos(theta_1 - theta_2))/(l1*l2^3*m2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));

      J0=[ddtheta1dm1 ddtheta2dm1; 
          ddtheta1dm2 ddtheta2dm2; 
          ddtheta1dl1 ddtheta2dl1; 
          ddtheta1dl2 ddtheta2dl2];

      % 计算在当前参数下实际输出和辨识输出差
      detlay(:,i) = y(:,i)-compute_output(:,i); 
      A = A + J0*J0';
      Q = Q - J0 * detlay(:,i); 
    end
    while true
              % 更新参数
        J_aim(t) = 0;
        v=-(A+(lambda/10)*I)^-1*Q;
        parameter(:,t) = parameter(:,t-1) + v ; 
        m1 = parameter(1,t); m2 = parameter(2,t); l1 = parameter(3,t); l2 = parameter(4,t);

        % 计算在当前参数下的目标函数
        for j = 1 : num

          theta_1 = x(1,j); theta_2 = x(2,j); p_1 = x(3,j); p_2 = x(4,j);

             dtheta_1(j)= (12*(2*l2*p_1 - 3*l1*p_2*cos(theta_1 - theta_2)))/(l1^2*l2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
   dtheta_2(j) = (12*(2*l1*m1*p_2 + 6*l1*m2*p_2 - 3*l2*m2*p_1*cos(theta_1 - theta_2)))/(l1*l2^2*m2*(8*m1 + 15*m2 - 9*m2*cos(2*theta_1 - 2*theta_2)));
   compute_output(:,j) = [dtheta_1(j);dtheta_2(j)];


             % 计算在当前参数下实际输出和辨识输出差
   detlay(:,j) = y(:,j)-compute_output(:,j); 
   % 根据当前参数计算目标函数值
   J_aim(t) = J_aim(t) +sum(power(detlay(:,j),2));
        end
        if J_aim(t) >= J_aim(t-1)
            lambda = 10*lambda;
        else 
            break;
        end
    end
    a = [abs(parameter(1,t)-parameter(1,t-1)) 
        abs(parameter(2,t)-parameter(2,t-1))
        abs(parameter(3,t)-parameter(3,t-1))
        abs(parameter(4,t)-parameter(4,t-1))];
    [m,p] =max(a);
    if m<e 
        break;
    else
        t=t+1;
        
        disp(J_aim(t));
    end
end
disp(parameter(:,t));





