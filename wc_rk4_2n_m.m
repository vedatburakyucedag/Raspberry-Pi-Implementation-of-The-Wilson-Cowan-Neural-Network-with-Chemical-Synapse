clc; clear; % 2 neuron wc standart deviation calculation

fns = 45; % fontsize

x_ele_sp1 = [rand rand];
y_ele_sp1 = [rand rand];
x_che_sp1 = [rand rand]; 
y_che_sp1 = [rand rand];

x_ele_sp2 = [rand rand]; 
y_ele_sp2 = [rand rand];
x_che_sp2 = [rand rand]; 
y_che_sp2 = [rand rand];

x_ele_sp3 = [rand rand];
y_ele_sp3 = [rand rand];
x_che_sp3 = [rand rand]; 
y_che_sp3 = [rand rand];

x_ele_sp4 = [rand rand];
y_ele_sp4 = [rand rand];
x_che_sp4 = [rand rand]; 
y_che_sp4 = [rand rand];

ns=50000;

ii = 1;
for kk=-5:0.01:5

    a_sp1 = 10; b_sp1 = 10; c_sp1 = 10; d_sp1 = 0;   Px_sp1 = 0.1; Py_sp1 = 0.1; Tx_sp1 = 0.8; Ty_sp1 = 0.8; Ix_sp1 = 0; Iy_sp1 = 0; u_sp1 = 0.5; % sp1
    a_sp2 = 10; b_sp2 = 10; c_sp2 = 10; d_sp2 = -10; Px_sp2 = 1;   Py_sp2 = 1;   Tx_sp2 = 0.3; Ty_sp2 = 0.3; Ix_sp2 = 0; Iy_sp2 = 0; u_sp2 = 0.3; % sp2
    a_sp3 = 2;  b_sp3 = 2;  c_sp3 = 5;  d_sp3 = -5;  Px_sp3 = -1;  Py_sp3 = 0.5; Tx_sp3 = 0.5; Ty_sp3 = 0.5; Ix_sp3 = 0; Iy_sp3 = 1; u_sp3 = 0.8; % sp3
    a_sp4 = 2;  b_sp4 = 2;  c_sp4 = 5;  d_sp4 = -5;  Px_sp4 = 1;   Py_sp4 = -3;  Tx_sp4 = 0.5; Ty_sp4 = 0.5; Ix_sp4 = 0; Iy_sp4 = 1; u_sp4 = 0.8; % sp4
    
    k = 10; Vs = -2; Teta_s = 0.28; c12 = 1 ;  t(1) = 0; wij_asenk = kk; h = 0.1; N = 500000;    
    
    [x_ele_sp1, y_ele_sp1, t] = function_wc_rk4_2n_m_ele(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_ele_sp1,y_ele_sp1,wij_asenk);
    [x_ele_sp2, y_ele_sp2, t] = function_wc_rk4_2n_m_ele(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_ele_sp2,y_ele_sp2,wij_asenk);
    [x_ele_sp3, y_ele_sp3, t] = function_wc_rk4_2n_m_ele(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_ele_sp3,y_ele_sp3,wij_asenk);
    [x_ele_sp4, y_ele_sp4, t] = function_wc_rk4_2n_m_ele(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_ele_sp4,y_ele_sp4,wij_asenk);
    [x_che_sp1, y_che_sp1, t] = function_wc_rk4_2n_m_che(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_che_sp1,y_che_sp1,wij_asenk,k,Vs,Teta_s,c12);
    [x_che_sp2, y_che_sp2, t] = function_wc_rk4_2n_m_che(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_che_sp2,y_che_sp2,wij_asenk,k,Vs,Teta_s,c12);
    [x_che_sp3, y_che_sp3, t] = function_wc_rk4_2n_m_che(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_che_sp3,y_che_sp3,wij_asenk,k,Vs,Teta_s,c12);
    [x_che_sp4, y_che_sp4, t] = function_wc_rk4_2n_m_che(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_che_sp4,y_che_sp4,wij_asenk,k,Vs,Teta_s,c12);

    R_ele(1,ii)=mean(std(x_ele_sp1(:,N-ns:N),1));
    R_ele(2,ii)=mean(std(x_ele_sp2(:,N-ns:N),1));
    R_ele(3,ii)=mean(std(x_ele_sp3(:,N-ns:N),1));
    R_ele(4,ii)=mean(std(x_ele_sp4(:,N-ns:N),1));
    R_che(1,ii)=mean(std(x_che_sp1(:,N-ns:N),1));
    R_che(2,ii)=mean(std(x_che_sp2(:,N-ns:N),1));
    R_che(3,ii)=mean(std(x_che_sp3(:,N-ns:N),1));
    R_che(4,ii)=mean(std(x_che_sp4(:,N-ns:N),1));

%     sttd1=sqrt((sum(x_ele_sp1(1,N-10000:N).^2)/10000 - (sum(x_ele_sp1(1,N-10000:N))/10000)^2)/(10000-1));
%     sttd2=sqrt((sum(x_che_sp1(2,N-10000:N).^2)/10000 - (sum(x_che_sp1(2,N-10000:N))/10000)^2)/(10000-1));
%     RR(1,ii)=mean([sttd1 sttd2]);

    wij(ii) = kk;
    ii = ii + 1;     

end

% Standard deviation plot

fig1 = figure('Position',get(0,'Screensize'));
plot(wij,R_ele(1,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig2 = figure('Position',get(0,'Screensize'));
plot(wij,R_ele(2,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig3 = figure('Position',get(0,'Screensize'));
plot(wij,R_ele(3,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig4 = figure('Position',get(0,'Screensize'));
plot(wij,R_ele(4,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig5 = figure('Position',get(0,'Screensize'));
plot(wij,R_che(1,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig6 = figure('Position',get(0,'Screensize'));
plot(wij,R_che(2,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig7 = figure('Position',get(0,'Screensize'));
plot(wij,R_che(3,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


fig8 = figure('Position',get(0,'Screensize'));
plot(wij,R_che(4,:),'LineStyle','-','Marker','o','Color','k','MarkerSize',25,'LineWidth',1)
ylabel('Standard Deviation (\sigma)')
xlabel('Snaptic Weight (g_{s})');
set(gca,'Fontsize',fns);
grid on


saveas(fig1, 'WC_ele_sp1_std.jpg');
saveas(fig2, 'WC_ele_sp2_std.jpg');
saveas(fig3, 'WC_ele_sp3_std.jpg');
saveas(fig4, 'WC_ele_sp4_std.jpg');
saveas(fig5, 'WC_che_sp1_std.jpg');
saveas(fig6, 'WC_che_sp2_std.jpg');
saveas(fig7, 'WC_che_sp3_std.jpg');
saveas(fig8, 'WC_che_sp4_std.jpg');
%%
clc; clear; % 2 neuron wilson cowan

%% asenkron
x_ele_sp1 = [-0.2230 0.2762];
y_ele_sp1 = [0.4162 -0.2942];
x_che_sp1 = [0.0520 -0.0985]; 
y_che_sp1 = [-0.6009 0.5841];

x_ele_sp2 = [-0.8144 0.8203]; 
y_ele_sp2 = [0.6778 0.1372];
x_che_sp2 = [0.85740 -0.8987]; 
y_che_sp2 = [0.3570 -0.0041];

x_ele_sp3 = [-0.8200 -0.4534];
y_ele_sp3 = [0.8555 -0.9568];
x_che_sp3 = [-0.9654 0.5173]; 
y_che_sp3 = [0.1229 -0.8311];

x_ele_sp4 = [0.9826 0.5382];
y_ele_sp4 = [0.3811 -0.7339];
x_che_sp4 = [0.9858 -0.2229]; 
y_che_sp4 = [-0.8523 0.9647];

%% senkron
x_ele_sp1 = [-0.6595 -0.6612];
y_ele_sp1 = [-0.4658 -0.4979];
x_che_sp1 = [-0.6358 -0.6483]; 
y_che_sp1 = [-0.3668 -0.4049];

x_ele_sp2 = [-0.8069 -0.8364]; 
y_ele_sp2 = [0.7066 0.5608];
x_che_sp2 = [0.3659 0.4629]; 
y_che_sp2 = [-0.8482 -0.8330];

x_ele_sp3 = [-0.5868 -0.6255];
y_ele_sp3 = [0.7732 0.7906];
x_che_sp3 = [0.5861 0.5547]; 
y_che_sp3 = [-0.1367 -0.0288];

x_ele_sp4 = [0.9705 0.9732];
y_ele_sp4 = [-0.8886 -0.8851];
x_che_sp4 = [0.1050 0.0400]; 
y_che_sp4 = [0.9780 0.9791];
% 
%%
a_sp1 = 10; b_sp1 = 10; c_sp1 = 10; d_sp1 = 0;   Px_sp1 = 0.1; Py_sp1 = 0.1; Tx_sp1 = 0.8; Ty_sp1 = 0.8; Ix_sp1 = 0; Iy_sp1 = 0; u_sp1 = 0.5; % sp1
a_sp2 = 10; b_sp2 = 10; c_sp2 = 10; d_sp2 = -10; Px_sp2 = 1;   Py_sp2 = 1;   Tx_sp2 = 0.3; Ty_sp2 = 0.3; Ix_sp2 = 0; Iy_sp2 = 0; u_sp2 = 0.3; % sp2
a_sp3 = 2;  b_sp3 = 2;  c_sp3 = 5;  d_sp3 = -5;  Px_sp3 = -1;  Py_sp3 = 0.5; Tx_sp3 = 0.5; Ty_sp3 = 0.5; Ix_sp3 = 0; Iy_sp3 = 1; u_sp3 = 0.8; % sp3
a_sp4 = 2;  b_sp4 = 2;  c_sp4 = 5;  d_sp4 = -5;  Px_sp4 = 1;   Py_sp4 = -3;  Tx_sp4 = 0.5; Ty_sp4 = 0.5; Ix_sp4 = 0; Iy_sp4 = 1; u_sp4 = 0.8; % sp4

k = 10; Vs = -2; Teta_s = 0.28; c12 = 1 ;  t(1) = 0; wij_asenk = -9; h = 0.01; N = 500000;  

[x_ele_sp1, y_ele_sp1, t] = function_wc_rk4_2n_m_ele(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_ele_sp1,y_ele_sp1,wij_asenk);
[x_ele_sp2, y_ele_sp2, t] = function_wc_rk4_2n_m_ele(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_ele_sp2,y_ele_sp2,wij_asenk);
[x_ele_sp3, y_ele_sp3, t] = function_wc_rk4_2n_m_ele(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_ele_sp3,y_ele_sp3,wij_asenk);
[x_ele_sp4, y_ele_sp4, t] = function_wc_rk4_2n_m_ele(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_ele_sp4,y_ele_sp4,wij_asenk);
[x_che_sp1, y_che_sp1, t] = function_wc_rk4_2n_m_che(a_sp1,b_sp1,c_sp1,d_sp1,Px_sp1,Py_sp1,Tx_sp1,Ty_sp1,Ix_sp1,Iy_sp1,u_sp1,h,N,t,x_che_sp1,y_che_sp1,wij_asenk,k,Vs,Teta_s,c12);
[x_che_sp2, y_che_sp2, t] = function_wc_rk4_2n_m_che(a_sp2,b_sp2,c_sp2,d_sp2,Px_sp2,Py_sp2,Tx_sp2,Ty_sp2,Ix_sp2,Iy_sp2,u_sp2,h,N,t,x_che_sp2,y_che_sp2,wij_asenk,k,Vs,Teta_s,c12);
[x_che_sp3, y_che_sp3, t] = function_wc_rk4_2n_m_che(a_sp3,b_sp3,c_sp3,d_sp3,Px_sp3,Py_sp3,Tx_sp3,Ty_sp3,Ix_sp3,Iy_sp3,u_sp3,h,N,t,x_che_sp3,y_che_sp3,wij_asenk,k,Vs,Teta_s,c12);
[x_che_sp4, y_che_sp4, t] = function_wc_rk4_2n_m_che(a_sp4,b_sp4,c_sp4,d_sp4,Px_sp4,Py_sp4,Tx_sp4,Ty_sp4,Ix_sp4,Iy_sp4,u_sp4,h,N,t,x_che_sp4,y_che_sp4,wij_asenk,k,Vs,Teta_s,c12);
%%
fig1 = figure('Position',get(0,'Screensize'));
plot(x_che_sp2(1,N-10000:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_che_sp2(2,N-10000:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp1 & E2 sp1')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp1','E2 sp1'},'Location','southwest');
% saveas(fig1, 'WC_ele_sp1_E1E2.jpg'); 

% save mat_wc_rk4_2n.mat x_ele_sp1 y_ele_sp1 x_ele_sp2 y_ele_sp2 x_ele_sp3 y_ele_sp3 x_ele_sp4 y_ele_sp4 x_che_sp1 y_che_sp1 x_che_sp2 y_che_sp2 x_che_sp3 y_che_sp3 x_che_sp4 y_che_sp4 

%%

fns = 65; % fontsize

ns1 = 1000;
ns2 = 1000;
ns3 = 1000;
ns4 = 1000;

% ns1 = N-1;
% ns2 = N-1;
% ns3 = N-1;
% ns4 = N-1;

fig1 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp1(1,N-ns1:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_ele_sp1(2,N-ns1:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp1 & E2 sp1')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp1','E2 sp1'},'Location','southwest');
saveas(fig1, 'WC_ele_sp1_E1E2.jpg'); 

fig2 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp2(1,N-ns2:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_ele_sp2(2,N-ns2:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp2 & E2 sp2')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp2','E2 sp2'},'Location','southwest');
saveas(fig2, 'WC_ele_sp2_E1E2.jpg'); 

fig3 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp3(1,N-ns3:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_ele_sp3(2,N-ns3:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp3 & E2 sp3')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp3','E2 sp3'},'Location','southwest');
saveas(fig3, 'WC_ele_sp3_E1E2.jpg'); 

fig4 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp4(1,N-ns4:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_ele_sp4(2,N-ns4:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp4 & E2 sp4')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp4','E2 sp4'},'Location','southwest');
saveas(fig4, 'WC_ele_sp4_E1E2.jpg'); 

fig5 = figure('Position',get(0,'Screensize'));
plot(x_che_sp1(1,N-ns1:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_che_sp1(2,N-ns1:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp1 & E2 sp1')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp1','E2 sp1'},'Location','southwest');
saveas(fig5, 'WC_che_sp1_E1E2.jpg'); 

fig6 = figure('Position',get(0,'Screensize'));
plot(x_che_sp2(1,N-ns2:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_che_sp2(2,N-ns2:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp2 & E2 sp2')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp2','E2 sp2'},'Location','southwest');
saveas(fig6, 'WC_che_sp2_E1E2.jpg'); 

fig7 = figure('Position',get(0,'Screensize'));
plot(x_che_sp3(1,N-ns3:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_che_sp3(2,N-ns3:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp3 & E2 sp3')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp3','E2 sp3'},'Location','southwest');
saveas(fig7, 'WC_che_sp3_E1E2.jpg'); 

fig8 = figure('Position',get(0,'Screensize'));
plot(x_che_sp4(1,N-ns4:N),'LineStyle','-','Marker','*','Color','r','MarkerSize',20,'LineWidth',8)
hold on
plot(x_che_sp4(2,N-ns4:N),'LineStyle','-','Marker','+','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp4 & E2 sp4')
xlabel('Number of Samples');
set(gca,'Fontsize',fns);
grid on
legend({'E1 sp4','E2 sp4'},'Location','southwest');
saveas(fig8, 'WC_che_sp4_E1E2.jpg'); 

%% --------------------------------------------------------------------------------------------------------------------------------
fns = 45; % fontsize

ns1 = N-2;
ns2 = N-2;
ns3 = N-2;
ns4 = N-2;

fig1 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp1(1,N-ns1:N),x_ele_sp1(2,N-ns1:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp1')
xlabel('E2 sp1');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig1, 'WC_ele_sp1_E1_E2.jpg'); 

fig2 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp2(1,N-ns2:N),x_ele_sp2(2,N-ns2:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp2')
xlabel('E2 sp2');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig2, 'WC_ele_sp2_E1_E2.jpg'); 

fig3 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp3(1,N-ns3:N),x_ele_sp3(2,N-ns3:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp3')
xlabel('E2 sp3');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig3, 'WC_ele_sp3_E1_E2.jpg'); 

fig4 = figure('Position',get(0,'Screensize'));
plot(x_ele_sp4(1,N-ns4:N),x_ele_sp4(2,N-ns4:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp4')
xlabel('E2 sp4');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig4, 'WC_ele_sp4_E1_E2.jpg'); 

fig5 = figure('Position',get(0,'Screensize'));
plot(x_che_sp1(1,N-ns1:N),x_che_sp1(2,N-ns1:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp1')
xlabel('E2 sp1');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig5, 'WC_che_sp1_E1_E2.jpg'); 

fig6 = figure('Position',get(0,'Screensize'));
plot(x_che_sp2(1,N-ns2:N),x_che_sp2(2,N-ns2:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp2')
xlabel('E2 sp2');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig6, 'WC_che_sp2_E1_E2.jpg'); 

fig7 = figure('Position',get(0,'Screensize'));
plot(x_che_sp3(1,N-ns3:N),x_che_sp3(2,N-ns3:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp3')
xlabel('E2 sp3');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig7, 'WC_che_sp3_E1_E2.jpg'); 

fig8 = figure('Position',get(0,'Screensize'));
plot(x_che_sp4(1,N-ns4:N),x_che_sp4(2,N-ns4:N),'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'LineWidth',8)
ylabel('E1 sp4')
xlabel('E2 sp4');
xlim([-2 2])
ylim([-1 1])
set(gca,'Fontsize',fns);
grid on
saveas(fig8, 'WC_che_sp4_E1_E2.jpg'); 