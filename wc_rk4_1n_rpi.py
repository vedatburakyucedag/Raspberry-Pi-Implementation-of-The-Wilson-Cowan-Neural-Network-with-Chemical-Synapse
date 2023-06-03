import math
import smbus

dac_address = [0x62 0x60 0x60 0x60]
mux_address = 0x70
mux_channel = [0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80]
bus = smbus.SMBus(1)
msg_zz = []

#x = -0.6098; y = -0.3053; a = 10; b = 10; c = 10; d = 0;   Px = 0.1; Py = 0.1; Tx = 0.8; Ty = 0.8; Ix = 0; Iy = 0; u = 0.5; # sp1
#x = 0.8700;  y = 0.5671;  a = 10; b = 10; c = 10; d = -10; Px = 1;   Py = 1;   Tx = 0.3; Ty = 0.3; Ix = 0; Iy = 0; u = 0.3; # sp2
#x = -0.7555; y = 0.8402;  a = 2;  b = 2;  c = 5;  d = -5;  Px = -1;  Py = 0.5; Tx = 0.5; Ty = 0.5; Ix = 0; Iy = 1; u = 0.8; # sp3
x = 0.9952;  y = -0.5925; a = 2;  b = 2;  c = 5;  d = -5;  Px = 1;   Py = -3;  Tx = 0.5; Ty = 0.5; Ix = 0; Iy = 1; u = 0.8; # sp4

h = 0.1

while 1:

        dx1 = (-x + math.tanh(u*(a*x - b*y + Px + Ix)))/Tx;
        dy1 = (-y + math.tanh(u*(c*x - d*y + Py + Iy)))/Ty;
        cx1 = x + (h/2) * dx1;
        cy1 = y + (h/2) * dy1;

        dx2 = (-cx1 + math.tanh(u*(a*cx1 - b*cy1 + Px + Ix)))/Tx;
        dy2 = (-cy1 + math.tanh(u*(c*cx1 - d*cy1 + Py + Iy)))/Ty;
        cx2 = x + (h/2) * dx2;
        cy2 = y + (h/2) * dy2;

        dx3 = (-cx2 + math.tanh(u*(a*cx2 - b*cy2 + Px + Ix)))/Tx;
        dy3 = (-cy2 + math.tanh(u*(c*cx2 - d*cy2 + Py + Iy)))/Ty;
        cx3 = x + h * dx3;
        cy3 = y + h * dy3;

        dx4 = (-cx3 + math.tanh(u*(a*cx3 - b*cy3 + Px + Ix)))/Tx;
        dy4 = (-cy3 + math.tanh(u*(c*cx3 - d*cy3 + Py + Iy)))/Ty;
        
        x = x + (h / 6) * (dx1 + 2 * dx2 + 2 * dx3 + dx4)
        y = y + (h / 6) * (dy1 + 2 * dy2 + 2 * dy3 + dy4)

        if x >= 0:

                zz = x * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[0])
                bus.write_i2c_block_data(dac_address[0], msg_f, msg_zz)

        else:

                zz = -x * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[1])
                bus.write_i2c_block_data(dac_address[1], msg_f, msg_zz)

        if y >= 0:

                zz = y * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[2])
                bus.write_i2c_block_data(dac_address[2], msg_f, msg_zz)

        else:

                zz = -y * 819.2
                msg_f = (int(zz) >> 8) & 0xff
                msg_s = (int(zz) & 0xff)
                msg_zz = [msg_s]
                bus.write_byte(mux_address, mux_channel[3])
                bus.write_i2c_block_data(dac_address[3], msg_f, msg_zz)


