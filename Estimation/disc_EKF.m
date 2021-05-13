function [xact,xhatp,Pp,T] = disc_EKF(x0,xhat0,P0,w,v,Q,R,r,q,t)
%This function simulates a system and implements a discrete extended kalman
%filter
%This function will need to be manually updated on a problem by problem
%basis
%x0 is the actual system initial conditions, nx1
%x_hat0 is the initial estimate of the system, nx1
%P0 is the initial estimate of the covariance matrix, nxn
%w is the ceoff matrix for the system noise, nx1
%v is the ceoff matrix for the measurement noise, m x 1
%q is the variance for the system noise
%r is the variance for the measurement noise
%t is the amount of time (in seconds to simulate the system)

%define constants used in problem
T=0.1; %seconds
num_steps=1:T:t;
N=[20 0];
E=[0 20];

%storage for important variables
xact=zeros(length(x0),length(num_steps));
xhatm=zeros(length(xhat0),length(num_steps));
xhatp=zeros(length(xhat0),length(num_steps));
[rp,cp]=size(P0);
Pm=zeros(rp,cp,length(num_steps));
Pp=zeros(rp,cp,length(num_steps));

%discrete EKF algorithm
for i=1:length(num_steps)
    %simulate actual system
    if i==1
        %initialization
        xact(:,i)=x0;
    else
        %noise
        W=normrnd(0,q)*w;
        V=normrnd(0,r)*v;
        %system matrices
        f=[1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
        h=[sqrt(((xact(1,i-1)-N(1))^2)+((xact(2,i-1)-E(1))^2));...
            sqrt(((xact(1,i-1)-N(2))^2)+((xact(2,i-1)-E(2))^2))];
        xact(:,i)=f*xact(:,i-1)+W;
        y=h+V;
    end
    if i==1
        %initialization
        xhatp(:,i)=xhat0;
        Pp(:,:,i)=P0;
    else
        %partial differential matrices
        F=f;
        L=eye(4,4);
        %covariance minus and xhat minus
        Pm(:,:,i)=F*Pp(:,:,i-1)*F'+L*Q*L';
        xhatm(:,i)=f*xhatp(:,i-1);
        h_xhatm=[sqrt(((xhatm(1,i)-N(1))^2)+((xhatm(2,i)-E(1))^2));...
            sqrt(((xhatm(1,i)-N(2))^2)+((xhatm(2,i)-E(2))^2))];
        h11=(xhatm(1,i)-N(1))/h_xhatm(1);
        h12=(xhatm(2,i)-E(1))/h_xhatm(1);
        h21=(xhatm(1,i)-N(2))/h_xhatm(2);
        h22=(xhatm(2,i)-E(2))/h_xhatm(2);
        H=[h11, h12 0 0;h21 h22 0 0];
        M=eye(2,2);
        %gain xhatp and Pp
        K=Pm(:,:,i)*H'/(H*Pm(:,:,i)*H'+M*R*M');
        xhatp(:,i)=xhatm(:,i)+K*(y-h_xhatm);
        Pp(:,:,i)=(eye(4,4)-K*H)*Pm(:,:,i);
    end
end

