function [x_che, y_che, phi_che, S_che, t] = function_wc_rk4_2n_m_che(a,b,c,d,Px,Py,Tx,Ty,Ix,Iy,u,h,N,t,x,y,wij,k1,k2,k3,alpha,beta,phi,dij,Vsyn,a0,To,Vshp,S)

nn = 2;

dx1 = zeros(1,nn);dy1 = zeros(1,nn);
dx2 = zeros(1,nn);dy2 = zeros(1,nn);
dx3 = zeros(1,nn);dy3 = zeros(1,nn);
dx4 = zeros(1,nn);dy4 = zeros(1,nn);
dphi1 = zeros(1,nn);
dphi2 = zeros(1,nn);
dphi3 = zeros(1,nn);
dphi4 = zeros(1,nn);
cphi1 = zeros(1,nn);
cphi2 = zeros(1,nn);
cphi3 = zeros(1,nn);

dS1 = zeros(1,nn);
dS2 = zeros(1,nn);
dS3 = zeros(1,nn);
dS4 = zeros(1,nn);
cS1 = zeros(1,nn);
cS2 = zeros(1,nn);
cS3 = zeros(1,nn);

cx1 = zeros(1,nn);cy1 = zeros(1,nn); 
cx2 = zeros(1,nn);cy2 = zeros(1,nn);
cx3 = zeros(1,nn);cy3 = zeros(1,nn);

coupling_phi = @(a3,a4,G) G*(a4-a3);

g = [0 wij ;
     wij 0 ];

G = [0 dij ;
     dij 0 ];

    for i =1:N

        t(i+1) = t(i) + h;

        for n = 1:nn

            if n == 2
                Isyn = g(n,n-1) * S(1,n-1) * (x(1,n) - Vsyn);
                dS1(1,n-1) = (a0 / (1 + exp(- x(1,n-1) / Vshp))) * (1 - S(1,n-1)) - (S(1,n-1) / To);
                cS1(1,n-1) = S(1,n-1) + (h / 2) * dS1(1,n-1);
                B = coupling_phi(phi(1,n), phi(1,n-1),G(n,n-1));
            else
                Isyn = g(n,n+1) * S(1,n+1) * (x(1,n) - Vsyn);
                dS1(1,n+1) = (a0 / (1 + exp(- x(1,n+1) / Vshp))) * (1 - S(1,n+1)) - (S(1,n+1) / To);
                cS1(1,n+1) = S(1,n+1) + (h / 2) * dS1(1,n+1);
                B = coupling_phi(phi(1,n), phi(1,n+1),G(n,n+1));
            end

            dx1(1,n) = (-x(1,n) + tanh(u*(a*x(1,n)-b*y(1,n)+Px+Ix-Isyn-k1*(alpha+3*beta*phi(1,n)^2)*x(1,n))))/Tx;
            dy1(1,n) = (-y(1,n) + tanh(u*(c*x(1,n)-d*y(1,n)+Py+Iy)))/Ty;
            dphi1(1,n) = k2*x(1,n) - k3*phi(1,n) + B;
            cx1(1,n) = x(1,n) + (h/2) * dx1(1,n);
            cy1(1,n) = y(1,n) + (h/2) * dy1(1,n);
            cphi1(1,n) = phi(1,n) + (h/2) * dphi1(1,n);

            if n == 2
                Isyn = g(n,n-1) * cS1(1,n-1) * (cx1(1,n) - Vsyn);
                dS2(1,n-1) = (a0 / (1 + exp(- cx1(1,n-1) / Vshp))) * (1 - cS1(1,n-1)) - (cS1(1,n-1) / To);
                cS2(1,n-1) = S(1,n-1) + (h / 2) * dS2(1,n-1);
                B = coupling_phi(cphi1(1,n), cphi1(1,n-1),G(n,n-1));
            else
                Isyn = g(n,n+1) * cS1(1,n+1) * (cx1(1,n) - Vsyn);
                dS2(1,n+1) = (a0 / (1 + exp(- cx1(1,n+1) / Vshp))) * (1 - cS1(1,n+1)) - (cS1(1,n+1) / To);
                cS2(1,n+1) = S(1,n+1) + (h / 2) * dS2(1,n+1);
                B = coupling_phi(cphi1(1,n), cphi1(1,n+1),G(n,n+1));
            end

            dx2(1,n) = (-cx1(1,n) + tanh(u*(a*cx1(1,n)-b*cy1(1,n)+Px+Ix-Isyn-k1*(alpha+3*beta*cphi1(1,n)^2)*cx1(1,n))))/Tx;
            dy2(1,n) = (-cy1(1,n) + tanh(u*(c*cx1(1,n)-d*cy1(1,n)+Py+Iy)))/Ty;
            dphi2(1,n) = k2*cx1(1,n) - k3*cphi1(1,n) + B;
            cx2(1,n) = x(1,n) + (h/2) * dx2(1,n);
            cy2(1,n) = y(1,n) + (h/2) * dy2(1,n);
            cphi2(1,n) = phi(1,n) + (h/2) * dphi2(1,n);

            if n == 2
                Isyn = g(n,n-1) * cS2(1,n-1) * (cx2(1,n) - Vsyn);
                dS3(1,n-1) = (a0 / (1 + exp(- cx2(1,n-1) / Vshp))) * (1 - cS2(1,n-1)) - (cS2(1,n-1) / To);
                cS3(1,n-1) = S(1,n-1) + h * dS3(1,n-1);
                B = coupling_phi(cphi2(1,n), cphi2(1,n-1),G(n,n-1));
            else
                Isyn = g(n,n+1) * cS2(1,n+1) * (cx2(1,n) - Vsyn);
                dS3(1,n+1) = (a0 / (1 + exp(- cx2(1,n+1) / Vshp))) * (1 - cS2(1,n+1)) - (cS2(1,n+1) / To);
                cS3(1,n+1) = S(1,n+1) + h * dS3(1,n+1);
                B = coupling_phi(cphi2(1,n), cphi2(1,n+1),G(n,n+1));
            end

            dx3(1,n) = (-cx2(1,n) + tanh(u*(a*cx2(1,n)-b*cy2(1,n)+Px+Ix-Isyn-k1*(alpha+3*beta*cphi2(1,n)^2)*cx2(1,n))))/Tx;
            dy3(1,n) = (-cy2(1,n) + tanh(u*(c*cx2(1,n)-d*cy2(1,n)+Py+Iy)))/Ty;
            dphi3(1,n) = k2*cx2(1,n) - k3*cphi2(1,n) + B;
            cx3(1,n) = x(1,n) + h * dx3(1,n);
            cy3(1,n) = y(1,n) + h * dy3(1,n);
            cphi3(1,n) = phi(1,n) + h * dphi3(1,n);

            if n == 2
                Isyn = g(n,n-1) * cS3(1,n-1) * (cx3(1,n) - Vsyn);
                dS4(1,n-1) = (a0 / (1 + exp(- cx3(1,n-1) / Vshp))) * (1 - cS3(1,n-1)) - (cS3(1,n-1) / To);
                B = coupling_phi(cphi3(1,n), cphi3(1,n-1),G(n,n-1));
            else
                Isyn = g(n,n+1) * cS3(1,n+1) * (cx3(1,n) - Vsyn);
                dS4(1,n+1) = (a0 / (1 + exp(- cx3(1,n+1) / Vshp))) * (1 - cS3(1,n+1)) - (cS3(1,n+1) / To);
                B = coupling_phi(cphi3(1,n), cphi3(1,n+1),G(n,n+1));
            end

            dx4(1,n) = (-cx3(1,n) + tanh(u*(a*cx3(1,n)-b*cy3(1,n)+Px+Ix-Isyn-k1*(alpha+3*beta*cphi3(1,n)^2)*cx3(1,n))))/Tx;
            dy4(1,n) = (-cy3(1,n) + tanh(u*(c*cx3(1,n)-d*cy3(1,n)+Py+Iy)))/Ty;
            dphi4(1,n) = k2*cx3(1,n) - k3*cphi3(1,n) + B;

            x(1,n) = x(1,n) + (h / 6) * (dx1(1,n) + 2 * dx2(1,n) + 2 * dx3(1,n) + dx4(1,n)); x_che(n,i+1) = x(1,n);
            y(1,n) = y(1,n) + (h / 6) * (dy1(1,n) + 2 * dy2(1,n) + 2 * dy3(1,n) + dy4(1,n)); y_che(n,i+1) = y(1,n);
            phi(1,n) = phi(1,n) + (h / 6) * (dphi1(1,n) + 2 * dphi2(1,n) + 2 * dphi3(1,n) + dphi4(1,n));  phi_che(n,i+1) = phi(1,n);
            S(1,n) = S(1,n) + (h / 6) * (dS1(1,n) + 2 * dS2(1,n) + 2 * dS3(1,n) + dS4(1,n)); S_che(n, i+1) = S(1,n);
            
        end
    end

end

