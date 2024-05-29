syms p_1 p_2 dtheta_1 dtheta_2  m1 m2 l1  l2  theta_1 theta_2 I1 I2 dx1 dx2 dy1 dy2 L g y1 y2 dp1 dp2 
I1=(1/12)*m1*l1^2;
I2 = (1/12)*m2*l2^2;
y1 = -(1/2)*l1*cos(theta_1);
dy1 = (1/2)*l1*sin(theta_1)*dtheta_1;

y2 =  -l1*cos(theta_1)-(1/2)*l2*cos(theta_2);
dy2 = l1*sin(theta_1)*dtheta_1 + (1/2)*l2*sin(theta_2)*dtheta_2;    

dx1 = (1/2)*l1*cos(theta_1)*dtheta_1;
dx2 = l1*cos(theta_1)*dtheta_1 + (1/2)*l2*cos(theta_2)*dtheta_2;

L = (1/2)*m1*(dx1^2+dy1^2) +(1/2)*m2*(dx2^2+dy2^2) +(1/2)*I1*(dtheta_1)^2 +(1/2)*I2*(dtheta_2)^2 -m1*g*y1-m1*g*y2;


[soldtheta_1, soldtheta_2] = solve([p_1==diff(L,'dtheta_1'),p_2 == diff(L,'dtheta_2')],[dtheta_1,dtheta_2]);

dp1=diff(L,'theta_1');
dp2=diff(L,'theta_2');

ddthe1dm1=diff(soldtheta_1,'m1');
ddthe1dm2=diff(soldtheta_1,'m2');
ddthe1dl1=diff(soldtheta_1,'l1');
ddthe1dl2=diff(soldtheta_1,'l2');

ddthe2dm1=diff(soldtheta_2,'m1');
ddthe2dm2=diff(soldtheta_2,'m2');
ddthe2dl1=diff(soldtheta_2,'l1');
ddthe2dl2=diff(soldtheta_2,'l2');

ddp1dm1=diff(dp1,'m1');
ddp1dm2=diff(dp1,'m2');
ddp1dl1=diff(dp1,'l1');
ddp1dl2=diff(dp1,'l2');

ddp2dm1=diff(dp2,'m1');
ddp2dm2=diff(dp2,'m2');
ddp2dl1=diff(dp2,'l1');
ddp2dl2=diff(dp2,'l2');
% 数学公式化处理
% dtheta1
stheta_1 = simplify(soldtheta_1);
ldtheta_1 = latex(stheta_1);
% dtheta2
sdtheta_2 = simplify(soldtheta_2);
ldtheta_2 = latex(sdtheta_2);
% dp1
sdp1 = simplify(dp1);
ldp1 = latex(sdp1);
% dp2
sdp2 = simplify(dp2);
ldp2 = latex(sdp2);

% ddthe1dm1、m2、l1、l2
sddthe1dm1 = simplify(ddthe1dm1);
lddthe1dm1 = latex(sddthe1dm1);
sddthe1dm2 = simplify(ddthe1dm2);
lddthe1dm2 = latex(sddthe1dm2);
sddthe1dl1 = simplify(ddthe1dl1);
lddthe1dl1 = latex(sddthe1dl1);
sddthe1dl2 = simplify(ddthe1dl2);
lddthe1dl2 = latex(sddthe1dl2);

% ddthe2dm1、m2、l1、l2
sddthe2dm1 = simplify(ddthe2dm1);
lddthe2dm1 = latex(sddthe2dm1);
sddthe2dm2 = simplify(ddthe2dm2);
lddthe2dm2 = latex(sddthe2dm2);
sddthe2dl1 = simplify(ddthe2dl1);
lddthe2dl1 = latex(sddthe2dl1);
sddthe2dl2 = simplify(ddthe2dl2);
lddthe2dl2 = latex(sddthe2dl2);

% ddp1dm1、m2、l1、l2
sddp1dm1 = simplify(ddp1dm1);
lddp1dm1 = latex(sddp1dm1);
sddp1dm2 = simplify(ddp1dm2);
lddp1dm2 = latex(sddp1dm2);
sddp1dl1 = simplify(ddp1dl1);
lddp1dl1 = latex(sddp1dl1);
sddp1dl2 = simplify(ddp1dl2);
lddp1dl2 = latex(sddp1dl2);

% ddp2dm1、m2、l1、l2
sddp2dm1 = simplify(ddp2dm1);
lddp2dm1 = latex(sddp2dm1);
sddp2dm2 = simplify(ddp2dm2);
lddp2dm2 = latex(sddp2dm2);
sddp2dl1 = simplify(ddp2dl1);
lddp2dl1 = latex(sddp2dl1);
sddp2dl2 = simplify(ddp2dl2);
lddp2dl2 = latex(sddp2dl2);