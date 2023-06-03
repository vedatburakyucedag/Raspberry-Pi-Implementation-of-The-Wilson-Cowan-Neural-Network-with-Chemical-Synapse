import math
import smbus

dac_address = [0x62, 0x62 0x60, 0x60]
mux_address = 0x70
mux_channel = [0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80]
bus = smbus.SMBus(1)
msg_zz = []

nn = 2

dx1 = []; dx2 = []; dx3 = []; dx4 = []
dy1 = []; dy2 = []; dy3 = []; dy4 = []
cx1 = []; cx2 = []; cx3 = []
cy1 = []; cy2 = []; cy3 = []

dx1 += [[0] * nn]; dx2 += [[0] * nn]; dx3 += [[0] * nn]; dx4 += [[0] * nn]
dy1 += [[0] * nn]; dy2 += [[0] * nn]; dy3 += [[0] * nn]; dy4 += [[0] * nn]
cx1 += [[0] * nn]; cx2 += [[0] * nn]; cx3 += [[0] * nn]
cy1 += [[0] * nn]; cy2 += [[0] * nn]; cy3 += [[0] * nn]

#a = 10; b = 10; c = 10; d = 0; Px = 0.1; Py = 0.1; Tx = 0.8; Ty = 0.8; Ix = 0; Iy = 0; u = 0.5; # sp1
#x = [-0.6639, 0.6460] # asenkron
#y = [-0.5741, 0.6160]
#x = [-0.6595, -0.6612]; # senkron
#y = [-0.4658, -0.4979];

#a = 10; b = 10; c = 10; d = -10; Px = 1; Py = 1; Tx = 0.3; Ty = 0.3; Ix = 0; Iy = 0; u = 0.3; # sp2
#x = [-0.9228, 0.8279] # asenkron
#y = [-0.4070, 0.7328]
#x = [-0.8069, -0.8364]; # senkron
#y = [0.7066, 0.5608];

#a = 2; b = 2; c = 5; d = -5; Px = -1; Py = 0.5; Tx = 0.5; Ty = 0.5; Ix = 0; Iy = 1; u = 0.8; # sp3
#x = [0.2415, -0.9565] # asenkron
#y = [-0.9789, 0.6342]
#x = [-0.5868, -0.6255]; # senkron
#y = [0.7732, 0.7906];

a = 2; b = 2; c = 5; d = -5; Px = 1; Py = -3; Tx = 0.5; Ty = 0.5; Ix = 0; Iy = 1; u = 0.8;  # sp4
#x = [0.5827, 0.9617] # asenkron
#y = [0.9443, -0.9064]
x = [0.9705, 0.9732]; # senkron
y = [-0.8886, -0.8851];

w2 = 0.01

g = [[0, w2], 
     [w2, 0]]
            
h = 0.1

while 1:

    for n in range(0, nn):

        if n == 1:
            A = g[n][n - 1] * (x[n - 1] - x[n])
        else:
            A = g[n][n + 1] * (x[n + 1] - x[n])

        dx1[0][n] = (-x[n] + math.tanh(u*(a*x[n] - b*y[n] + Px + Ix + A)))/Tx;
        dy1[0][n] = (-y[n] + math.tanh(u*(c*x[n] - d*y[n] + Py + Iy)))/Ty;
        cx1[0][n] = x[n] + (h/2) * dx1[0][n];
        cy1[0][n] = y[n] + (h/2) * dy1[0][n];

        if n == 1:
            A = g[n][n - 1] * (cx1[0][n - 1] - cx1[0][n])
        else:
            A = g[n][n + 1] * (cx1[0][n + 1] - cx1[0][n])

        dx2[0][n] = (-cx1[0][n] + math.tanh(u*(a*cx1[0][n] - b*cy1[0][n] + Px + Ix + A)))/Tx;
        dy2[0][n] = (-cy1[0][n] + math.tanh(u*(c*cx1[0][n] - d*cy1[0][n] + Py + Iy)))/Ty;
        cx2[0][n] = x[n] + (h/2) * dx2[0][n];
        cy2[0][n] = y[n] + (h/2) * dy2[0][n];
        
        if n == 1:
            A = g[n][n - 1] * (cx2[0][n - 1] - cx2[0][n])
        else:
            A = g[n][n + 1] * (cx2[0][n + 1] - cx2[0][n])

        dx3[0][n] = (-cx2[0][n] + math.tanh(u*(a*cx2[0][n] - b*cy2[0][n] + Px + Ix + A)))/Tx;
        dy3[0][n] = (-cy2[0][n] + math.tanh(u*(c*cx2[0][n] - d*cy2[0][n] + Py + Iy)))/Ty;
        cx3[0][n] = x[n] + h * dx3[0][n];
        cy3[0][n] = y[n] + h * dy3[0][n];

        if n == 1:
            A = g[n][n - 1] * (cx3[0][n - 1] - cx3[0][n])
        else:
            A = g[n][n + 1] * (cx3[0][n + 1] - cx3[0][n])

        dx4[0][n] = (-cx3[0][n] + math.tanh(u*(a*cx3[0][n] - b*cy3[0][n] + Px + Ix + A)))/Tx;
        dy4[0][n] = (-cy3[0][n] + math.tanh(u*(c*cx3[0][n] - d*cy3[0][n] + Py + Iy)))/Ty;

        x[n] = x[n] + (h / 6) * (dx1[0][n] + 2 * dx2[0][n] + 2 * dx3[0][n] + dx4[0][n])
        y[n] = y[n] + (h / 6) * (dy1[0][n] + 2 * dy2[0][n] + 2 * dy3[0][n] + dy4[0][n])
        
        if x[n] >= 0:

                zz = x[n] * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[n])
                bus.write_i2c_block_data(dac_address[n], msg_f, msg_zz)

        else:

                zz = -x[n] * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[n+2])
                
                bus.write_i2c_block_data(dac_address[n+2], msg_f, msg_zz)
