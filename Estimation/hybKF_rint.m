function [x_act,x_hat_p,P_p] = hybKF_rint(x0,x_hat0,P0,t,T)
%This function implements the hybrid Kalman filter algorithm on a system
%Due to the complex nature of the problems associated with nonlinear
%systems, this function will most likely need to be updated on a
%per-problem basis
%This function uses rectangular integration methods
%x0 is the actual system initial conditions, nx1
%x_hat0 is the initial estimate of the system, nx1
%P0 is the initial estimate of the covariance matrix, nxn
%w is the ceoff matrix for the system noise, nx1
%v is the ceoff matrix for the measurement noise, m x 1
%q is the variance for the system noise
%r is the variance for the measurement noise
%T is the integration time step
%t is the total time to simulate the system

%define number of time steps to integrate
num_steps=0:T:0.5;  %each measurement occurs every 0.5 seconds

%storage for important variables
x_act=zeros(length(x0),length(1:0.5:t));
x_hat_m=zeros(length(x_hat0),length(1:0.5:t));
x_hat_p=zeros(length(x_hat0),length(1:0.5:t));
[rp,cp]=size(P0);
P_m=zeros(rp,cp,length(1:0.5:t));
P_p=zeros(rp,cp,length(1:0.5:t));
K=zeros(rp,1,length(1:0.5:t));

%define constants needed for problem dynamics
rho0=2;
g=32.2;
k=20000;
R=10000;
Q=zeros(3,3);
M_h=100000;
a=100000;

%hybrid kalman filter algorithm
for i=1:length(1:0.5:t)
   %calculate random variables
    W=normrnd(0,sqrt(Q(1,1)),length(x_hat0));
    V=normrnd(0,sqrt(R),1);
    %Now we simulate the actual system and implement the filter
    %actual system
    if i==1
        %initialization
        x_act(:,i)=x0;
    else
        %otherwise
        x_act_d=x_act(:,i-1);
        for j1=1:length(num_steps)
            x_act_d=x_act_d(:)+T*[x_act_d(2)+W(1)... 
        rho0*exp(-x_act_d(1)/k)*(x_act_d(2)^2)*(x_act_d(3)/2)-g+W(2)...
        W(3)]';
        end
        x_act(:,i)=x_act_d;
    end
    %measurement
    y=sqrt(M_h^2+(x_act(1,i)-a)^2)+V;
    %implement the filter
    if i==1
        %initialization
        x_hat_p(:,i)=x_hat0;
        P_p(:,:,1)=P0;
    else
        %integrate state estimate and covariance matrix from k-1(+) to k(-)
        x_hat_d=x_hat_p(:,i-1);   %set dummy variable equal to xhatk-1(+)
        P_d=P_p(:,:,i-1);   %set dummy equal to Pk-1(+)
        %integrate to xhatk(-) and Pk(-)
        for j2=1:length(num_steps)
            %calculate differential matrices
            A21=(-rho0*exp(-x_hat_d(1)/k)*(x_hat_d(2)^2)*x_hat_d(3))/(2*k);
            A22=rho0*exp(-x_hat_d(1)/k)*x_hat_d(2)*x_hat_d(3);
            A23=(rho0*exp(-x_hat_d(1)/k)*(x_hat_d(2)^2))/(2);
            A=[0 1 0;A21 A22 A23;0 0 0];
            L=eye(3,3);
            C1=(x_hat_d(1)-a)/sqrt(M_h^2+(x_hat_d(1)-a)^2);
            C=[C1 0 0];
            M=1;
            % calculate Q_tild and R_tild
            Q_tild=L*Q*L';
            R_tild=M*R*M';
            %
            x_hat_d=x_hat_d+T*([x_hat_d(2)...
                rho0*exp(-x_hat_d(1)/k)*(x_hat_d(2)^2)*x_hat_d(3)/(2)-g...
                0]');
            P_d=P_d+T*(A*P_d+P_d*A'+Q_tild);
        end
        x_hat_m(:,i)=x_hat_d;   %set xhatk(-)
        P_m(:,:,i)=P_d; %set Pk(-)
        %calculate x_hatk(+), Pk(+) and Kk
        K(:,:,i)=P_m(:,:,i)*C'*(C*P_m(:,:,i)*C'+R_tild)^-1;
        x_hat_p(:,i)=x_hat_m(:,i)+K(:,i)*(y-sqrt(M_h^2+(x_hat_m(1,i)-a)^2));
        P_p(:,:,i)=(eye(rp,cp)-K(:,:,i)*C)*P_m(:,:,i)...
            *(eye(rp,cp)-K(:,:,i)*C)'+K(:,:,i)*R_tild*K(:,:,i)';
    end
end
end

