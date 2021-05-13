function [x,xhatp,Pp] = UKF(x0,xhatp0,Pp0,T,t)
%This function implements the unscented Kalman filter algorithm on a system
%due to the complex nature of problems, this algorithm will need to be
%updated on a problem by problem basis

%calculate number of time steps
num_steps=1:0.5:t;
update_step=0:T:0.5;

%create storage
sigpoints=zeros(length(xhatp0),2*length(xhatp0));
x=zeros(length(x0),length(num_steps));
xhatm=zeros(length(xhatp0),length(num_steps));
xhat=zeros(length(xhatp0),2*length(xhatp0),length(num_steps));
xhatp=zeros(length(xhatp0),length(num_steps));
y_temp=zeros(1,2*length(x0),length(num_steps));
yhat=zeros(1,length(num_steps));
y=zeros(1,length(num_steps));
Pm=zeros(length(xhatp0),length(xhatp0),length(num_steps));
Pp=zeros(length(xhatp0),length(xhatp0),length(num_steps));
K=zeros(length(xhatp0),1,length(num_steps));

%Problem specific values
%constants
p0=2;
g=32.2;
kk=20000;   %do not get confused with variable k used in the loop
R=10000;
Q=zeros(3,3);
M=100000;
a=100000;
n=length(xhatp0);

%initialize system
W=normrnd(0,0,length(xhatp0),1);
x_d=x0;
for i=1:length(update_step)
    f=[x_d(2);...
        p0*exp(-x_d(1)/kk)*(x_d(2)^2)*x_d(3)/2-g;...
        0];
    x_d=x_d+T*(f+W);
end
x(:,1)=x_d;

%UKF algorithm
for k=1:length(num_steps)
    %calculate random variables
    W=normrnd(0,0,length(xhatp0),1);
    V=normrnd(0,sqrt(R),1);
    %actual system
    x_d=x(:,k);
    for j=1:length(update_step)
        f=[x_d(2);...
            p0*exp(-x_d(1)/kk)*(x_d(2)^2)*x_d(3)/2-g;...
            0];
        x_d=x_d+T*(f+W);
    end
    x(:,k+1)=x_d;
    % measurement
    y(:,k)=(M^2+(x(1,k)-a)^2)^(1/2)+V;
    %time update equations
    %sigma points
    if k==1
        %sigma points
        [root, p] = chol(n*Pp0);
        sigpoints=xhatp0+root';
        sigpoints=[sigpoints xhatp0-root'];
    else
        %sigma points
        [root,p] = chol(n*Pp(:,:,k-1));
        sigpoints=xhatp(:,k-1)+root';
        sigpoints=[sigpoints xhatp(:,k-1)-root'];
    end
   %transform sigma points
   for i=1:2*n
       xhat_d=sigpoints(:,i);
       for j=1:length(update_step)
           fxhat=[xhat_d(2);...
                p0*exp(-xhat_d(1)/kk)*(xhat_d(2)^2)*xhat_d(3)/2-g;...
                0];
            xhat_d=xhat_d+T*fxhat;
       end
       xhat(:,i,k)=xhat_d;
   end
    %combine vectors to get a priori state estimate
    apsum=0;
    for i=1:2*n
        apsum=apsum+(1/(2*n))*xhat(:,i,k);
    end
    xhatm(:,k)=apsum;
    %a priori error covariance
    sum2=zeros(n,n);
    for i=1:2*n
        sum2=sum2+(1/(2*n))*((xhat(:,i,k)-xhatm(:,k))*(xhat(:,i,k)-xhatm(:,k))');
    end 
    Pm(:,:,k)=sum2+Q;
    %measurement update equations
    %use original sigma points at time k
    %calculate predicted measurement
    for i=1:2*n
        y_temp(:,i,k)=(M^2+(xhat(1,i,k)-a)^2)^(1/2);
    end
    sum3=0;
    for i=1:2*n
        sum3=sum3+(1/(2*n))*y_temp(:,i,k);
    end
    yhat(:,k)=sum3;
    %calculate Py
    sum4=0;
    for i=1:2*n
        sum4=sum4+(1/(2*n))*(y_temp(:,i,k)-yhat(:,k))*...
            (y_temp(:,i,k)-yhat(:,k))';
    end
    Py=sum4+R;
    %calculate Pxy
    sum5=0;
    for i=1:2*n
        sum5=sum5+(1/(2*n))*((xhat(:,i,k)-xhatm(:,k))*...
            (y_temp(:,i,k)-yhat(:,k))');
    end
    Pxy=sum5;
    %calculate gain
    K(:,:,k)=Pxy*(Py)^-1;
    %calculate xhatp and Pp
    xhatp(:,k)=xhatm(:,k)+K(:,:,k)*(y(:,k)-yhat(:,k));
    Pp(:,:,k)=Pm(:,:,k)-K(:,:,k)*Py*K(:,:,k)';
end
end