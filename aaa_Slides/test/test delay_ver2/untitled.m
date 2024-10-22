%set measure real time > signal 'a'


trig.ioObj = io64;
%initialize the inpoutx64 system driver
status = io64(trig.ioObj);
if(status == 0)
    disp("EMG triggers ready")
end
trig.address = portAddress;

io64(trig.ioObj, trig.address, 1);
io64(trig.ioObj, trig.address, 0);

while 1
    if 'a' <-1
        toc
        break
    end
end
