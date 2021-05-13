%Nathan Lutes
%MADAM DNN test - 2 link planar robot arm
%addpath('Z:\Desktop\Graduate Research\Spr 21\MultWeightUpdate')
addpath('S:\Desktop\Graduate Research\Spr 21\MultWeightUpdate')
clear; clc; close all
global tauhist y r

%constants and initial conditions
states = 4; in = states; out = 4; L = 5; tauhist = [];
gbarV = 0; gbarW1 = 0; gbarW2 = 0; gbarW3 = 0; gbarWout = 0;
%eta = 0.01; etaP = 3*eta; sigPr = 10; B = 0.999;
eta = 0.001; etaP = 100*eta; sigPr = 1250; B = 0.999;
%etaP = 100*eta;

x0 = [0.5; 0.5; 0; 0];
xhat0 = [0.5,0.5,0,0]';
xti = x0;
V = normrnd(0,3.6/(sqrt(L)),[in+1,L]);
W1 = normrnd(0,3.6/(sqrt(L)),[L+1,L]);
W2 = normrnd(0,3.6/(sqrt(L)),[L+1,L]);
W3 = normrnd(0,3.6/(sqrt(L)),[L+1,L]);
Wout = normrnd(0,3.6/(sqrt(L)),[L+1,out]);

%encode weights
weights.V = V;
weights.W1 = W1;
weights.W2 = W2;
weights.W3 = W3;
weights.Wout = Wout;

%simulation
dt = 0.001;
t0 = 0;
tf = 10;
t1 = t0:dt:tf;
x1 = [x0; xhat0];
%storage
x1hist = zeros(length(x1),length(t0:dt:tf));

%algorithm
for i = 1:length(t1)
    %save current state
    x1hist(:,i) = x1;
    %calculate neural network and gradient
    y = NNOut(weights,x1(1:4));
    %calculate gradient
    [xdot] = RobInvDynControl(t1(i),x1);
    
    for j = 1:1
        %weight update
        %madam algorithm
        g = getGrad(weights,x1(1:4));
        gV = g.gV;
        gW1 = g.gW1;
        gW2 = g.gW2;
        gW3 = g.gW3;
        gWout = g.gWout;
        gbarV = sqrt((1-B)*gV.^2 + B*gbarV.^2);
        gbarW1 = sqrt((1-B)*gW1.^2 + B*gbarW1.^2);
        gbarW2 = sqrt((1-B)*gW2.^2 + B*gbarW2.^2);
        gbarW3 = sqrt((1-B)*gW3.^2 + B*gbarW3.^2);
        gbarWout = sqrt((1-B)*gWout.^2 + B*gbarWout.^2);
        gcV = clamp(gV./gbarV, -(etaP/eta),etaP/eta);
        gcW1 = clamp(gW1./gbarW1, -(etaP/eta),etaP/eta);
        gcW2 = clamp(gW2./gbarW2, -(etaP/eta),etaP/eta);
        gcW3 = clamp(gW3./gbarW3, -(etaP/eta),etaP/eta);
        gcWout = clamp(gWout./gbarWout, -(etaP/eta),etaP/eta);
        Vi = V .* exp(-eta*sign(V).*gcV);
        W1i = W1 .* exp(-eta*sign(W1).*gcW1);
        W2i = W2 .* exp(-eta*sign(W2).*gcW2);
        W3i = W3 .* exp(-eta*sign(W3).*gcW3);
        Wouti = Wout .* exp(-eta*sign(Wout).*gcWout);
        V = clamp(Vi,-sigPr,sigPr);
        W1 = clamp(W1i,-sigPr,sigPr);
        W2 = clamp(W2i,-sigPr,sigPr);
        W3 = clamp(W3i,-sigPr,sigPr);
        Wout = clamp(Wouti,-sigPr,sigPr);
        %encode weights
        weights.V = V;
        weights.W1 = W1;
        weights.W2 = W2;
        weights.W3 = W3;
        weights.Wout = Wout;
    end
    
    x1 = x1 + dt*xdot;
end

%calculate desired trajectory
w = 0.5;
xd = [sin(w*t1); cos(w*t1)];

%plot desired and actual trajectories
Fsize = 15;
Pic_Width=7;
Pic_Height=7;
% figure
% hold on;
% plot(t,x1,'g','LineWidth',3);
% plot(t,x2,':m','LineWidth',3);
% hold off
% xlabel('\bf{Time(s)}','FontWeight','bold','FontSize',20,'interpreter','latex')
% ylabel('\bf{$x_1$ (rad)}','FontWeight','bold','FontSize',20,'interpreter','latex')
% ax = gca;ax.LineWidth = 3; ax.FontSize =Fsize; box on; ax.FontWeight='bold';
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 Pic_Width Pic_Height]);
% grid;set(gca,'MinorGridLineStyle','-');set(gca,'GridLineStyle','-.');box on;

figure
plot(t1,x1hist(1,:),'LineWidth',3)
hold on
plot(t1,x1hist(2,:),'LineWidth',3)
plot(t1,xd(1,:),'--','LineWidth',3)
plot(t1,xd(2,:),'--','LineWidth',3)
hold off
legend('q1','q2','qd1','qd2')
xlabel('\bf{Time(s)}','FontWeight','bold','FontSize',20,'interpreter','latex')
ylabel('\bf{Angle (rad)}','FontWeight','bold','FontSize',20,'interpreter','latex')
ax = gca;ax.LineWidth = 3; ax.FontSize =Fsize; box on; ax.FontWeight='bold';
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 Pic_Width Pic_Height]);
grid;set(gca,'MinorGridLineStyle','-');set(gca,'GridLineStyle','-.');box on;

%plot error 
figure
e = x1hist(1:2,:) - xd;
plot(t1,e(1,:),'LineWidth',3);
hold on
plot(t1,e(2,:),'LineWidth',3);
hold off
legend('e1','e2')
xlabel('\bf{Time(s)}','FontWeight','bold','FontSize',20,'interpreter','latex')
ylabel('\bf{Angle (rad)}','FontWeight','bold','FontSize',20,'interpreter','latex')
ax = gca;ax.LineWidth = 3; ax.FontSize =Fsize; box on; ax.FontWeight='bold';
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 Pic_Width Pic_Height]);
grid;set(gca,'MinorGridLineStyle','-');set(gca,'GridLineStyle','-.');box on;

%plot tau history
figure
plot(t1,tauhist(1,:),'LineWidth',3)
hold on
plot(t1,tauhist(2,:),'LineWidth',3)
hold off
legend('tau1','tau2')
xlabel('\bf{Time(s)}','FontWeight','bold','FontSize',20,'interpreter','latex')
ylabel('\bf{Torque}','FontWeight','bold','FontSize',20,'interpreter','latex')
ax = gca;ax.LineWidth = 3; ax.FontSize =Fsize; box on; ax.FontWeight='bold';
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 Pic_Width Pic_Height]);
grid;set(gca,'MinorGridLineStyle','-');set(gca,'GridLineStyle','-.');box on;
print('Torq','-depsc','-r300')


% %plot event count
% figure
% plot(t1,evtcnt)
% xlabel('Time')
% ylabel('Number of Events vs. Time')
% title('Event count history')
