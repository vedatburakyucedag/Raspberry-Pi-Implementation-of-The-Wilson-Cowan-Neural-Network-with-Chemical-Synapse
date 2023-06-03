clc; clear; % 2 neuron fhn


k1 = 0; k2=0.9; k3=0.5; alpha=0.4; beta=0.02; dij = 0.2; % k1 = 0; alan etkisi yok, k1 sifirdan farkli, alan etkisi var.

phi_sp1 = [1 -1]; 
phi_sp2 = [-1 1];  
phi_sp3 = [1 -1];  
phi_sp4 = [1 -1];  

a_sp1 = 10; b_sp1 = 10; c_sp1 = 10; d_sp1 = 0; Px_sp1 = 0.1; Py_sp1 = 0.1; Tx_sp1 = 0.8; Ty_sp1 = 0.8; Ix_sp1 = 0; Iy_sp1 = 0; u_sp1 = 0.5; % sp1
x_ele_sp1 = [-0.6639 0.6460]; % asenkron
y_ele_sp1 = [-0.5741 0.6160];
x_che_sp1 = [-0.0373 -0.0166]; 
y_che_sp1 = [-0.6916 0.6844];

% x_ele_sp1 = [-0.659590384822278 -0.661234588002797]; % senkron
% y_ele_sp1 = [-0.465811718110040 -0.497993447202711];
% x_che_sp1 = [-0.635856661641105 -0.648376282934035]; 
% y_che_sp1 = [-0.366823774169340 -0.404978648913966];

a_sp2 = 10; b_sp2 = 10; c_sp2 = 10; d_sp2 = -10; Px_sp2 = 1; Py_sp2 = 1; Tx_sp2 = 0.3; Ty_sp2 = 0.3; Ix_sp2 = 0; Iy_sp2 = 0; u_sp2 = 0.3; % sp2
x_ele_sp2 = [-0.9228 0.8279]; % asenkron
y_ele_sp2 = [-0.4070 0.7328];
x_che_sp2 = [0.8836 -0.9199]; 
y_che_sp2 = [0.4621 -0.1591];

% x_ele_sp2 = [-0.806919437611357 -0.836444470699158]; % senkron
% y_ele_sp2 = [0.706688345106675 0.560838225602429];
% x_che_sp2 = [0.365962140928225 0.462925183827205]; 
% y_che_sp2 = [-0.848278872041117 -0.833012975534743];

a_sp3 = 2; b_sp3 = 2; c_sp3 = 5; d_sp3 = -5; Px_sp3 = -1; Py_sp3 = 0.5; Tx_sp3 = 0.5; Ty_sp3 = 0.5; Ix_sp3 = 0; Iy_sp3 = 1; u_sp3 = 0.8; % sp3
x_ele_sp3 = [0.2415 -0.9565]; % asenkron
y_ele_sp3 = [-0.9789 0.6342];
x_che_sp3 = [-0.5430 -0.8223]; 
y_che_sp3 = [-0.9490 0.8640];

% x_ele_sp3 = [-0.586814148870581 -0.625572521933029]; % senkron
% y_ele_sp3 = [0.773299059158284 0.790689182897620];
% x_che_sp3 = [0.586134389437051 0.554724056415338]; 
% y_che_sp3 = [-0.136711823080494 -0.0288241858803498];

a_sp4 = 2; b_sp4 = 2; c_sp4 = 5; d_sp4 = -5; Px_sp4 = 1; Py_sp4 = -3; Tx_sp4 = 0.5; Ty_sp4 = 0.5; Ix_sp4 = 0; Iy_sp4 = 1; u_sp4 = 0.8;  % sp4
x_ele_sp4 = [0.5827 0.9617]; % asenkron
y_ele_sp4 = [0.9443 -0.9064];
x_che_sp4 = [0.858657178480275 0.826597261672139]; 
y_che_sp4 = [-0.884545305528627 0.824698682257132];

% x_ele_sp4 = [0.970520464065232 0.973223290847922]; % senkron
% y_ele_sp4 = [-0.888683875189944 -0.885114084354817];
% x_che_sp4 = [0.105093193207749 0.0400627036270377]; 
% y_che_sp4 = [0.978042614667581 0.979145099778186];

Vsyn = 0; a0 = 2; To = 5/6; Vshp = 0.05; t(1) = 0; wij_asenk = -0.2;
S_sp1 = [0.189359003160682 0.606027057486384];
S_sp2 = [0.606107382034013 0.163666026975214];  
S_sp3 = [0.591336334456183 0.0374639722270593]; 
S_sp4 = [0.602155728249883 0.624998580000941];

h = 0.1; N = 300;

[x_ele_sp1, y_ele_sp1, phi_ele_sp1, t] = function_wc_rk4_2n_m_ele(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_ele_sp1,y_ele_sp1,wij_asenk,k1,k2,k3,alpha,beta,phi_sp1,dij);
[x_ele_sp2, y_ele_sp2, phi_ele_sp2, t] = function_wc_rk4_2n_m_ele(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_ele_sp2,y_ele_sp2,wij_asenk,k1,k2,k3,alpha,beta,phi_sp2,dij);
[x_ele_sp3, y_ele_sp3, phi_ele_sp3, t] = function_wc_rk4_2n_m_ele(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_ele_sp3,y_ele_sp3,wij_asenk,k1,k2,k3,alpha,beta,phi_sp3,dij);
[x_ele_sp4, y_ele_sp4, phi_ele_sp4, t] = function_wc_rk4_2n_m_ele(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_ele_sp4,y_ele_sp4,wij_asenk,k1,k2,k3,alpha,beta,phi_sp4,dij);
[x_che_sp1, y_che_sp1, phi_che_sp1, S_che_sp1, t] = function_wc_rk4_2n_m_che(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_che_sp1,y_che_sp1,wij_asenk,k1,k2,k3,alpha,beta,phi_sp1,dij,Vsyn,a0,To,Vshp,S_sp1);
[x_che_sp2, y_che_sp2, phi_che_sp2, S_che_sp2, t] = function_wc_rk4_2n_m_che(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_che_sp2,y_che_sp2,wij_asenk,k1,k2,k3,alpha,beta,phi_sp2,dij,Vsyn,a0,To,Vshp,S_sp2);
[x_che_sp3, y_che_sp3, phi_che_sp3, S_che_sp3, t] = function_wc_rk4_2n_m_che(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_che_sp3,y_che_sp3,wij_asenk,k1,k2,k3,alpha,beta,phi_sp3,dij,Vsyn,a0,To,Vshp,S_sp3);
[x_che_sp4, y_che_sp4, phi_che_sp4, S_che_sp4, t] = function_wc_rk4_2n_m_che(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_che_sp4,y_che_sp4,wij_asenk,k1,k2,k3,alpha,beta,phi_sp4,dij,Vsyn,a0,To,Vshp,S_sp4);

% A = 0;
% 
% for i=1:length(x_ele_sp1)
% 
%     dd(i) = sqrt((x_ele_sp3(1,i) - x_ele_sp3(2,i))^2 + (y_ele_sp3(1,i) - y_ele_sp3(2,i))^2);
% 
%     A = A + dd;
% 
% end
% 
% mean(dd)
% 
% A/length(x_ele_sp1)
% 
% dd = sqrt((x_ele_sp1(1,1000) - x_ele_sp1(2,1000))^2 + (y_ele_sp1(1,1000) - y_ele_sp1(2,1000))^2);

% save mat_wc_rk4_2n.mat x_ele_sp1 y_ele_sp1 x_ele_sp2 y_ele_sp2 x_ele_sp3 y_ele_sp3 x_ele_sp4 y_ele_sp4 x_che_sp1 y_che_sp1 x_che_sp2 y_che_sp2 x_che_sp3 y_che_sp3 x_che_sp4 y_che_sp4 

% figure(1); clf(1);
% plot(x_ele_sp1(1,:),'k','LineWidth',1)
% hold on
% plot(x_ele_sp1(2,:),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(2); clf(2);
% plot(x_ele_sp2(1,:),'k','LineWidth',1)
% hold on
% plot(x_ele_sp2(2,:),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(3); clf(3);
% plot(x_ele_sp3(1,:),'k','LineWidth',1)
% hold on
% plot(x_ele_sp3(2,:),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(4); clf(4);
% plot(x_ele_sp4(1,:),'k','LineWidth',1)
% hold on
% plot(x_ele_sp4(2,:),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on

figure(5); clf(5);
plot(x_che_sp1(1,:),'k','LineWidth',1)
hold on
plot(x_che_sp1(2,:),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(6); clf(6);
plot(x_che_sp2(1,:),'k','LineWidth',1)
hold on
plot(x_che_sp2(2,:),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(7); clf(7);
plot(x_che_sp3(1,:),'k','LineWidth',1)
hold on
plot(x_che_sp3(2,:),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(8); clf(8);
plot(x_che_sp4(1,:),'k','LineWidth',1)
hold on
plot(x_che_sp4(2,:),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

%%

% figure(1); clf(1);
% plot(x_ele_sp1(1,250000:N),'k','LineWidth',1)
% hold on
% plot(x_ele_sp1(2,250000:N),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(2); clf(2);
% plot(x_ele_sp2(1,250000:N),'k','LineWidth',1)
% hold on
% plot(x_ele_sp2(2,250000:N),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(3); clf(3);
% plot(x_ele_sp3(1,250000:N),'k','LineWidth',1)
% hold on
% plot(x_ele_sp3(2,250000:N),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on
% 
% figure(4); clf(4);
% plot(x_ele_sp4(1,250000:N),'k','LineWidth',1)
% hold on
% plot(x_ele_sp4(2,250000:N),'r','LineWidth',1)
% xlabel('time')
% ylabel('x')
% set(gca,'Fontsize',10)
% grid on

figure(5); clf(5);
plot(x_che_sp1(1,250000:N),'k','LineWidth',1)
hold on
plot(x_che_sp1(2,250000:N),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(6); clf(6);
plot(x_che_sp2(1,250000:N),'k','LineWidth',1)
hold on
plot(x_che_sp2(2,250000:N),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(7); clf(7);
plot(x_che_sp3(1,250000:N),'k','LineWidth',1)
hold on
plot(x_che_sp3(2,250000:N),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on

figure(8); clf(8);
plot(x_che_sp4(1,250000:N),'k','LineWidth',1)
hold on
plot(x_che_sp4(2,250000:N),'r','LineWidth',1)
xlabel('time')
ylabel('x')
set(gca,'Fontsize',10)
grid on