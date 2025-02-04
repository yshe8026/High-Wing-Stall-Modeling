clear all; close all;
tic
u = udp('127.0.0.1',49000,'LocalPortMode','manual','LocalPort',49747);
% u.DatagramTerminateMode = 'on';

fopen(u);
readasync(u);

fwrite(u,['RPOS', 0x00, '5', 0x00]);

T = 2; %record time in seconds
T_finish = now()+T/86400;
n = 0;
while now()<=T_finish
    if u.BytesAvailable >0
        if u.BytesAvailable>69
            disp(u.BytesAvailable);
        end
        A = fread(u,69);
        toc
        tic
        n = n+1;
    end
%     size(A)
end
disp(n)
disp(u.BytesAvailable)

fclose(u);