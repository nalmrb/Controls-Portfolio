%Nathan Lutes
%ctrlsspecialassign2
%12/1/2018

%Part a
A=[-0.49 0.5e-4 -4.79 0; 0 -0.148e-1 -13.877 -32.2; 1 -0.19e-3 -0.836 0; 1 0 0 0];
B=[-8.743; -1.096; -0.01115; 0];
C=zeros(4,4);
D=zeros(4,1);
Q=[14.79 0 0 0; 0 2.7e-10 1.182e-6 0; 0 1.182e-6 0.00514 0; 0 0 0 10.41];
R=51.02+9.101e-7;
N=[0; 1.57e-8; 6.85e-5; 0];
T=linspace(0,50,10000);

%calculate open loop pole locations
sysO=ss(A, B, C, D);
Eopen=pole(sysO);
Eopen

%calculate gains and closed loop poles
[Kcl, Scl, Ecl]=lqr(A,B,Q,R,N);
Acl=A-B*Kcl;
syscl=ss(Acl,B,C,D);
Kcl
Ecl

%part B response plot
%Initial condition 1
[yopen,t,xopen]=initial(sysO,[0;0;5/57.3;0],T);
[ycl,tcl,xcl]=initial(syscl,[0;0;5/57.3;0],T);

%q
plot(T,xopen(:,1))
hold on
plot(T,xcl(:,1))
hold off
legend('q open','q closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('q Response vs. Time for Condition 1')
%v
figure()
plot(T,xopen(:,2))
hold on
plot(T,xcl(:,2))
hold off
legend('v open','v closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('v Response vs. Time for Condition 1')
%a
figure()
plot(T,xopen(:,3))
hold on
plot(T,xcl(:,3))
hold off
legend('a open','a closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('a Response vs. Time for Condition 1')
%theta
figure()
plot(T,xopen(:,4))
hold on
plot(T,xcl(:,4))
hold off
legend('theta open','theta closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('theta Response vs. Time for Condition 1')

%part B response plot
%Initial conditions 2

[yopen2,t2,xopen2]=initial(sysO,[0;0;5/57.3;5/57.3],T);
[ycl2,tcl2,xcl2]=initial(syscl,[0;0;5/57.3;5/57.3],T);

%q
figure()
plot(T,xopen2(:,1))
hold on
plot(T,xcl2(:,1))
hold off
legend('q open','q closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('q Response vs. Time for Condition 2')
%v
figure()
plot(T,xopen2(:,2))
hold on
plot(T,xcl2(:,2))
hold off
legend('v open','v closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('v Response vs. Time for Condition 2')
%a
figure()
plot(T,xopen2(:,3))
hold on
plot(T,xcl2(:,3))
hold off
legend('a open','a closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('a Response vs. Time for Condition 2')
%theta
figure()
plot(T,xopen2(:,4))
hold on
plot(T,xcl2(:,4))
hold off
legend('theta open','theta closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('theta Response vs. Time for Condition 2')

%Part C(a)
Qc=[41.62 0 0 0; 0 2.7e-10 1.182e-6 0; 0 1.182e-6 0.00514 0; 0 0 0 10.41];

%calculate open loop pole locations
sysOc=ss(A, B, C, D);
Eopenc=pole(sysOc);
Eopenc

%calculate gains and closed loop poles
[Kclc, Sclc, Eclc]=lqr(A,B,Qc,R,N);
Aclc=A-B*Kclc;
sysclc=ss(Aclc,B,C,D);
Kclc
Eclc

%Part C(B)
%initial conditions 1
[yopenc,tc,xopenc]=initial(sysOc,[0;0;5/57.3;0],T);
[yclc,tclc,xclc]=initial(sysclc,[0;0;5/57.3;0],T);

%q
figure()
plot(T,xopenc(:,1))
hold on
plot(T,xclc(:,1))
hold off
legend('q open','q closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('q Response Part C vs. Time for Condition 1')
%v
figure()
plot(T,xopenc(:,2))
hold on
plot(T,xclc(:,2))
hold off
legend('v open','v closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('v Response Part Cvs. Time for Condition 1')
%a
figure()
plot(T,xopenc(:,3))
hold on
plot(T,xclc(:,3))
hold off
legend('a open','a closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('a Response Part C vs. Time for Condition 1')
%theta
figure()
plot(T,xopenc(:,4))
hold on
plot(T,xclc(:,4))
hold off
legend('theta open','theta closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('theta Response Part C vs. Time for Condition 1')

%part B response plot
%Initial conditions 2

[yopen2c,t2c,xopen2c]=initial(sysOc,[0;0;5/57.3;5/57.3],T);
[ycl2c,tcl2c,xcl2c]=initial(sysclc,[0;0;5/57.3;5/57.3],T);

%q
figure()
plot(T,xopen2c(:,1))
hold on
plot(T,xcl2c(:,1))
hold off
legend('q open','q closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('q Response Part C vs. Time for Condition 2')
%v
figure()
plot(T,xopen2c(:,2))
hold on
plot(T,xcl2c(:,2))
hold off
legend('v open','v closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('v Response Part Cvs. Time for Condition 2')
%a
figure()
plot(T,xopen2c(:,3))
hold on
plot(T,xcl2c(:,3))
hold off
legend('a open','a closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('a Response Part C vs. Time for Condition 2')
%theta
figure()
plot(T,xopen2c(:,4))
hold on
plot(T,xcl2c(:,4))
hold off
legend('theta open','theta closed')
xlabel('Time (seconds)')
ylabel('System Response')
title('theta Response Part C vs. Time for Condition 2')