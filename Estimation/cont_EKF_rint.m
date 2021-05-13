function [x_act,x_hat,P,K] = cont_EKF_rint(x0,x_hat0,P0,t,T)
%This function calculates the continuous kalman filter for a given system
%This function uses rectangular integration to simulate the filter and the
%system
%initial conditions
%x0 is the actual system initial conditions, nx1
%x_hat0 is the initial estimate of the system, nx1
%P0 is the initial estimate of the covariance matrix, nxn
%w is the ceoff matrix for the system noise, nx1
%v is the ceoff matrix for the measurement noise, m x 1
%q is the variance for the system noise
%r is the variance for the measurement noise
%t is the amount of time (in seconds to simulate the system)
%T is the time step with which to conduct the rectangular integration

%define number of time steps
num_steps=1:T:t;

%create storage of important variables
[rP,cP]=size(P0);
P=zeros(rP,cP,length(num_steps));
x_act=zeros(length(x0),length(num_steps));
x_hat=zeros(length(x_hat0),length(num_steps));
K=zeros(length(x0),1,length(num_steps));

%define constants needed for problem dynamics
rho0=2;
g=32.2;
k=20000;
R=10000;
Q=zeros(3,3);
M_h=100000;
a=100000;

%continuous Kalman filter algorithm
%initialize
x_act(:,1)=x0;
x_hat(:,1)=x_hat0;
P(:,:,1)=P0;
C1=(x_hat(1,1)-a)/sqrt(M_h^2+(x_hat(1,1)-a)^2);
C=[C1 0 0];
M=1;
R_tild=M*R*M';
K(:,:,1)=P(:,:,1)*C'/R_tild;
for i=1:length(num_steps)
    %calculate random variables
    W=normrnd(0,sqrt(Q(1,1)),length(x_hat0));
    V=normrnd(0,sqrt(R),1);
    %simulate real system and calculate measurement
    %due to the nonlinear nature of the problem, the system dynamics will 
    %need to be hardcoded for each particular problem
    %propogate actual system
    x_act(:,i+1)=x_act(:,i)+T*[x_act(2,i)+W(1)... 
        rho0*exp(-x_act(1,i)/k)*(x_act(2,i)^2)*(x_act(3,i)/2)-g+W(2)...
        W(3)]';
    %next measurement
    y=sqrt(M_h^2+(x_act(1,i)-a)^2)+V;
    %The filter algorithm
    %calculate differential matrices
    A21=(-rho0*exp(-x_hat(1,i)/k)*(x_hat(2,i)^2)*x_hat(3,i))/(2*k);
    A22=rho0*exp(-x_hat(1,i)/k)*x_hat(2,i)*x_hat(3,i);
    A23=(rho0*exp(-x_hat(1,i)/k)*(x_hat(2,i)^2))/(2);
    A=[0 1 0;A21 A22 A23;0 0 0];
    L=eye(3,3);
    C1=(x_hat(1,i)-a)/sqrt(M_h^2+(x_hat(1,i)-a)^2);
    C=[C1 0 0];
    M=1;
    % calculate Q_tild and R_tild
    Q_tild=L*Q*L';
    R_tild=M*R*M';
    %calculate estimate
    %calculate gain
    K(:,:,i+1)=P(:,:,i)*C'/R_tild;
    x_hat(:,i+1)=x_hat(:,i)+T*([x_hat(2,i)...
        (rho0*exp(-x_hat(1,i)/k)*(x_hat(2,i)^2)*(x_hat(3,i)/2))-g...
        0]'+K(:,:,i)*(y-sqrt(M_h^2+(x_hat(1,i)-a)^2)));
    %calculate covariance matrix
    P(:,:,i+1)=P(:,:,i)+T*(A*P(:,:,i)+P(:,:,i)*A'+Q_tild-P(:,:,i)...
        *C'*inv(R_tild)*C*P(:,:,i));
end
end
