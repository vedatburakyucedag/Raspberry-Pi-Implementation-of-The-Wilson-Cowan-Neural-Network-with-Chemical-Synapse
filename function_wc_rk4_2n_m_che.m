function [x_asenk,y_asenk,t] = function_wc_rk4_2n_m_che(a,b,c,d,Px,Py,Tx,Ty,Ix,Iy,u,h,N,t,x,y,wij,k,Vs,Teta_s,c12)

nn = 2;

dx1 = zeros(1,nn);dy1 = zeros(1,nn);
dx2 = zeros(1,nn);dy2 = zeros(1,nn);
dx3 = zeros(1,nn);dy3 = zeros(1,nn);
dx4 = zeros(1,nn);dy4 = zeros(1,nn);

cx1 = zeros(1,nn);cy1 = zeros(1,nn); 
cx2 = zeros(1,nn);cy2 = zeros(1,nn);
cx3 = zeros(1,nn);cy3 = zeros(1,nn);

coupling = @(a1,a2,g) g * (a1 - Vs) * c12 * (1 / (1 + exp(- k * (a2 - Teta_s))));
B=0;
g = [0 wij ;
     wij 0 ];

    for i =1:N

        t(i+1) = t(i) + h;

        for n = 1:nn

            if n == 2
                A = coupling(x(1,n), x(1,n - 1), g(n,n - 1));
%                 B = coupling(y(1,n), y(1,n - 1), g(n,n - 1));
            else
                A = coupling(x(1,n), x(1,n + 1), g(n,n + 1));
%                 B = coupling(y(1,n), y(1,n + 1), g(n,n + 1));
            end

            dx1(1,n) = (-x(1,n) + tanh(u*(a*x(1,n)-b*y(1,n)+Px+Ix+A)))/Tx;
            dy1(1,n) = (-y(1,n) + tanh(u*(c*x(1,n)-d*y(1,n)+Py+Iy+B)))/Ty;
            cx1(1,n) = x(1,n) + (h/2) * dx1(1,n);
            cy1(1,n) = y(1,n) + (h/2) * dy1(1,n);

            if n == 2
                A = coupling(cx1(1,n), cx1(1,n - 1), g(n,n - 1));
%                 B = coupling(cy1(1,n), cy1(1,n - 1), g(n,n - 1));
            else
                A = coupling(cx1(1,n), cx1(1,n + 1), g(n,n + 1));
%                 B = coupling(cy1(1,n), cy1(1,n + 1), g(n,n + 1));
            end

            dx2(1,n) = (-cx1(1,n) + tanh(u*(a*cx1(1,n)-b*cy1(1,n)+Px+Ix+A)))/Tx;
            dy2(1,n) = (-cy1(1,n) + tanh(u*(c*cx1(1,n)-d*cy1(1,n)+Py+Iy+B)))/Ty;
            cx2(1,n) = x(1,n) + (h/2) * dx2(1,n);
            cy2(1,n) = y(1,n) + (h/2) * dy2(1,n);

            if n == 2
                A = coupling(cx2(1,n), cx2(1,n - 1), g(n,n - 1));
%                 B = coupling(cy2(1,n), cy2(1,n - 1), g(n,n - 1));
            else
                A = coupling(cx2(1,n), cx2(1,n + 1), g(n,n + 1));
%                 B = coupling(cy2(1,n), cy2(1,n + 1), g(n,n + 1));
            end

            dx3(1,n) = (-cx2(1,n) + tanh(u*(a*cx2(1,n)-b*cy2(1,n)+Px+Ix+A)))/Tx;
            dy3(1,n) = (-cy2(1,n) + tanh(u*(c*cx2(1,n)-d*cy2(1,n)+Py+Iy+B)))/Ty;
            cx3(1,n) = x(1,n) + h * dx3(1,n);
            cy3(1,n) = y(1,n) + h * dy3(1,n);

            if n == 2
                A = coupling(cx3(1,n), cx3(1,n - 1), g(n,n - 1));
%                 B = coupling(cy3(1,n), cy3(1,n - 1), g(n,n - 1));
            else
                A = coupling(cx3(1,n), cx3(1,n + 1), g(n,n + 1));
%                 B = coupling(cy3(1,n), cy3(1,n + 1), g(n,n + 1));
            end

            dx4(1,n) = (-cx3(1,n) + tanh(u*(a*cx3(1,n)-b*cy3(1,n)+Px+Ix+A)))/Tx;
            dy4(1,n) = (-cy3(1,n) + tanh(u*(c*cx3(1,n)-d*cy3(1,n)+Py+Iy+B)))/Ty;

            x(1,n) = x(1,n) + (h / 6) * (dx1(1,n) + 2 * dx2(1,n) + 2 * dx3(1,n) + dx4(1,n)); x_asenk(n,i+1) = x(1,n);
            y(1,n) = y(1,n) + (h / 6) * (dy1(1,n) + 2 * dy2(1,n) + 2 * dy3(1,n) + dy4(1,n)); y_asenk(n,i+1) = y(1,n);

        end
    end

end

