function [x,u_xhatp,xhatp,P] = particle_filter(x0,t)
%This function implements a particle filter given a set of initial
%conditions. Due to the complexity of most problems, this algorithm will
%need to be updated on a per-problem basis
%note this function is currently setup for a scalar problem

%calculate number of time steps
num_steps=1:t;

%choose particles
%set number of particles, choose 500 for this problem
N=500;

%create storage
%note that the storage was setup for a scalar problem and that it will need
%to be modified for a vector problem by implementing 3D matrices
x=zeros(length(x0),length(num_steps));
y=zeros(1,length(num_steps));
xhatm=zeros(N,length(num_steps));
xhatp=zeros(N,length(num_steps));
q=zeros(N,length(num_steps));
qs=zeros(N,length(num_steps));
u_xhatp=zeros(1,length(num_steps));
P=zeros(1,length(num_steps));

%useful variables
%m=number of measurement equation dimensions
m=1;
R=1;
Q=10;

%Initial pdf of X
u_xhatp0=0;
var_xhatp0=2;

%generate initial particles
xhatp0=normrnd(u_xhatp0,sqrt(var_xhatp0),N,1);

%particle filter algorithm
for k=1:t
    %random variables
    W=normrnd(0,sqrt(Q),1,1);
    V=normrnd(0,sqrt(R),1,1);
    %actual system simulation
    if k==1
        %initialize system
        uk=8*cos(1.2*(num_steps(k)-1));
        x(:,1)=(0.5*x0+(25*x0)/(1+x0^2)+uk+W);
    else
        uk=8*cos(1.2*(num_steps(k)-1));
        x(:,k)=(0.5*x(:,k-1)+(25*x(:,k-1))/(1+x(:,k-1)^2)+uk+W);
    end
    %measurement
    y(:,k)=((x(:,k)^2)/20)+V;
    %time propogation
    if k==1
        for i=1:N
            %generate noise for every particle
            W=normrnd(0,sqrt(Q),1,1);
            xhatm(i,k)=(0.5*xhatp0(i)+(25*xhatp0(i))/(1+xhatp0(i)^2)+uk+W);
        end
    else
        for i=1:N
            %generate noise for every particle
            W=normrnd(0,sqrt(Q),1,1);
            xhatm(i,k)=(0.5*xhatp(i,k-1)+...
                (25*xhatp(i,k-1))/(1+xhatp(i,k-1)^2)+uk+W);
        end
    end
    %compute relative likelihood qi of each particle
    for i=1:N
        h=((xhatm(i,k)^2)/20);
        q1=((2*pi)^(m/2)*abs(R)^(1/2))^-1;
        q2=exp(0.5*(-(y(:,k)-h)'*(R^-1)*(y(:,k)-h)));
        q(i,k)=q1*q2;
    end
    %scale relatively likelihoods
    qsum=0;
    for i=1:N
        qsum=qsum+q(i,k);
    end
    qs(:,k)=q(:,k)/qsum;
    %generate a posteriori particles on the basis of the relative
    %likelihoods (resampling step)
    for i=1:N
        qsum2=0;
        r=rand; %uniform random number on interval (0,1)
        for j=1:N
            qsum2=qsum2+qs(j,k);
            if qsum2<abs(r)
                continue
            else
                break
            end
        end
        xhatp(i,k)=xhatm(j,k);
    end
    %compute desired statistical measurements
    %compute mean of particles
    u_xhatp(k)=mean(xhatp(:,k));
    %calculate sample covariance
    cov_sum=0;
    for i=1:N
        cov=(xhatp(i,k)-u_xhatp(k))*(xhatp(i,k)-u_xhatp(k))';
        cov_sum=cov_sum+cov;
    end
    P(:,k)=cov_sum/N;
end

