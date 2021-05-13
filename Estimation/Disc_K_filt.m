function [P_plus,P_minus,x_act,x_hat_plus,x_hat_minus,y,K] = Disc_K_filt(x,P,x_hat,F,H,w,v,q,r,t)
%This function computes the discrete kalman filter for a state space system
%Note that this function does not include control
%noise is assumed to have mean zero and a specified variance
%It is intended to be generalized for many types of problems.
%x is the initial state vector, it is an nx1 vector
%P is the initial cov matrix, it is an nxn matrix
%x_hat is the initial state estimate, nxn
%F is the state dynamics matrix, n x n
%H is the measurement coefficient matrix, m x n
%w is the noise coeff matrix associated with the state model, n x 1 matrix
%v is the noise coeff matrix associated with the measurement, m x 1 matrix
%q is the variance of w
%r is the variance of v
%t is the number of time steps to run

% initialize storage of important variables
[r,c]=size(P);
P_minus=zeros(r,c,t);
P_plus=zeros(r,c,t);
x_act=zeros(length(x),t);
x_hat_minus=zeros(length(x_hat),t);
x_hat_plus=zeros(length(x_hat),t);
[rH,cH]=size(H);
y = zeros(rH,t);
K=zeros(r,rH,t);
%add initial values
P_plus(:,:,1)=P;
x_hat_plus(:,1)=x_hat;
x_act(:,1)=x;
%create Q and R matrices
Q=eye(length(w),length(w))*q;
R=eye(length(v),length(v))*r;

% Kalman filter algorithm
for i=1:t
    if i==1
        x_act(:,i)=x;
    else
       %calculate random variables
       W=normrnd(0,q,size(w)).*w;
       V=normrnd(0,r,size(v)).*v;
       %state propogation 
       x_act(:,i)=F*x_act(:,i-1)+W; %calculate actual state
       y(:,i)=H*x_act(:,i) + V;     %calculate measurement
    end
    if i==1
        x_hat_plus(:,i)=x_hat;
        P_plus(:,:,i)=P;
    else
       x_hat_minus(:,i)=F*x_hat_plus(:,i-1);    %propogate state using model
       P_minus(:,:,i)=F*P_plus(:,:,i-1)*F'+Q;   %propogate covariance
       K(:,:,i)=P_minus(:,:,i)*H'/((H*P_minus(:,:,i)*H' + R));...
           %calculate gain matrix
       x_hat_plus(:,i)=x_hat_minus(:,i) + K(:,:,i)*(y(:,i) - H*x_hat_minus(:,i));...
           %calculate new estimate using measurement
       P_plus(:,:,i)=(eye(r,r)-K(:,:,i)*H)*P_minus(:,:,i)...
           *(eye(r,r)-K(:,:,i)*H)'+K(:,:,i)*R*K(:,:,i)';
    end
end
end

