
function dxdt=ode_func(t,x)
global A B q1 q2 q1dot1 q2dot1 delm1 delm2 q1dot2 q2dot2 t_old; 
global delm11 delm21 kp1 kp2 kd1 kd2 m1 m2 a1 a2 g time1 j m1hat m2hat  V Vdot ;

global m1predict;
global m2predict;
m1predict = m1+ (delm1);
m2predict = m2+ (delm2);
time1(j)=t;

Mp=M(m1predict,m2predict,q2);
Cp=C(m1predict,m2predict,q2,q1dot1,q2dot1);
Np=N(m1predict,m2predict,q1,q2);
     
M1=M(m1,m2,q2);  
C1=C(m1,m2,q2,q1dot1,q2dot1);
N1=N(m1,m2,q1,q2);

kp=[kp1 0;0 kp2];
kd=[kd1 0;0 kd2];

q=[x(1);x(2)];
dq=[x(3);x(4)];

tau = Cp*dq + Np +Mp*([-sin(t);cos(t)] -kp*(q-[sin(t);-cos(t)])-kd*(dq-[cos(t);sin(t)]));
W =vpa([q1dot2*(a1^2)+g*a1*cos(q1),q1dot2*(a1^2+a2^2+2*a1*a2*cos(q2))+q2dot2*(a2^2 + a1*a2*cos(q2))-a1*a2*sin(q2)*q2dot1*(2*q1dot1 + q2dot1)+g*a1*cos(q1)+g*a2*cos(q1+q2);0,q1dot2*(a2^2 + a1*a2*cos(q2))+q2dot2*(a2^2)+a1*a2*sin(q2)*(q1dot1^2)+g*a2*cos(q1 + q2)]);  

      
dxdt = zeros(6,1);
dxdt(1) = x(3);
dxdt(2) = x(4);
%k=[-sin(t);cos(t)] -kp*(q-[sin(t);-cos(t)])-kd*(dq-[cos(t);sin(t)])+ Mp\W*[x(5);x(6)];
k = M1\(tau-C1*[dxdt(1);dxdt(2)]-N1) ;
dxdt(3) = k(1,1);
dxdt(4) = k(2,1);

display(t);
W =vpa([q1dot2*(a1^2)+g*a1*cos(q1),q1dot2*(a1^2+a2^2+2*a1*a2*cos(q2))+q2dot2*(a2^2 + a1*a2*cos(q2))-a1*a2*sin(q2)*q2dot1*(2*q1dot1 + q2dot1)+g*a1*cos(q1)+g*a2*cos(q1+q2);0,q1dot2*(a2^2 + a1*a2*cos(q2))+q2dot2*(a2^2)+a1*a2*sin(q2)*(q1dot1^2)+g*a2*cos(q1 + q2)]);  
Mp=M(m1predict,m2predict,q2);

qd1= sin(t);
qd2 = -cos(t);
Dqd11=cos(t);
Dqd21=sin(t);
Dqd12 =-sin(t);
Dqd22 = cos(t);
Q=[1 0 0 0;0 1 0 0; 0 0 1 0;0 0 0 1];
A1 = [0 0 1 1;0 0 1 1; -kp1 0 -kd1 0; 0 -kp2 0 -kd2];
P1 = lyap(A',Q);


k2= -2*transpose(W)*transpose(inv(Mp))*(B')*P1*([x(1)-qd1;x(2)-qd2;x(3)-Dqd11;x(4)-Dqd21]);

time1(j)=t;

dxdt(5) = k2(1,1);
dxdt(6) = k2(2,1);
    
delm1 =x(5);
delm2 =x(6);   
q1 =x(1); 
q2 = x(2);
q1dot1 =x(3);
q2dot1=x(4);
q1dot2=dxdt(3);
q2dot2=dxdt(4);


m1hat(j) =m1predict;
m2hat(j)=m2predict;
j=j+1;
delX =[x(1)-qd1;x(2)-qd2;x(3)-Dqd11;x(4)-Dqd21];
deltheta= [delm1;delm2];
V(j)=delX'*P1*delX + 0.5*(deltheta')*deltheta;
Vdot(j) = delX'*(P1*A +A'*P1)*delX +2*deltheta'*W'*inv(Mp')*B'*P1*delX + deltheta'*k2;


end