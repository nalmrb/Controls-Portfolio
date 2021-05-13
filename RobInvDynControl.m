function [qdot] = RobInvDynControl(t,q)
%This function implements an event triggered neural adaptive controller for
%use in NN project

%define event triggered state
global tau tauhist y r 

%constants
states = 4; in = states;
K2 = [30*eye(2),0*eye(2);zeros(2), 30*eye(2)]; 
a1 = 1; a2 = 1; m1 = 1; m2 = 2.3; g = 9.8;

%get input measurements
x = q(1:in); xhat = q(in+1:length(q));

%calculate desired trajectory
w = 0.5;
xd = [sin(w*t); cos(w*t)];
xddot = [w*cos(w*t); -w*sin(w*t)];
F = [0 0 1 0; 0 0 0 1; -w^2 0 0 0; 0 -w^2 0 0]*[xd; xddot];

%Robot system
M = [(m1+m2)*a1^2 + m2*a2^2 + 2*m2*a1*a2*cos(x(2))...
    m2*a2^2 + m2*a1*a2*cos(x(2));...
    m2*a2^2 + m2*a1*a2*cos(x(2)) m2*a2^2];
N = [-m2*a1*a2*(2*x(3)*x(4) + x(4)^2)*sin(x(2)) + (m1+m2)*g*a1*cos(x(1)) + m2*g*a2*cos(x(1) + x(2));...
    m2*a1*a2*(x(3)^2)*sin(x(2)) + m2*g*a2*cos(x(1) + x(2))];
B = [zeros(2);inv(M)];
Bs = [-inv(M);inv(M)];
Us = [1;1];
Baug = [B,Bs];

%tracking error
%r = -([xd;xddot]-[x(1); x(2); x(3); x(4)]);
r = x-xhat;

%periodic
%calculate control input
%y = 0;
K = [7.5*eye(2), 0*eye(2);0*eye(2),7.5*eye(2)];
tau = Baug\(F - y - K*(x-[xd;xddot]) + Bs*Us);
%Baug\(F - K*(xhat-[xd;xddot]) - K2*(xti-xhat) - y + Bs*Us);

%record control input
tauhist = [tauhist tau];
% %record event history
% evtcnt = [evtcnt cnt];

%calculate state and weight update laws
xdot = [x(3); x(4); -M\N] + Baug*tau - Bs*Us;% + B*d;
xhatdot = y + Baug*tau - Bs*Us + K2*(x-xhat);

qdot = [xdot; xhatdot];

end