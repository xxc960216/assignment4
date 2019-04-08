clear all

runTime = 1; %given in seconds
timecuts = 1000;
dt =runTime/timecuts;

R1=1;
C1=0.25;
R2=2;
L1=0.2;
R3=10;
a=100;
R4=0.1;
Ro=1000;
Cn=0.00001;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo
   
   %     V1             V2  I1    V3  I3     V4  In            Vo  
G = [     1,             0,  0,    0,  0,     0,  0,            0; ... V1
      -1/R1, (1/R2 + 1/R1), -1,    0,  0,     0,  0,            0; ... V2
          0,             1,  0,   -1,  0,     0,  0,            0; ... I1
          0,             0, -1, 1/R3,  0,     0, -1,            0; ... V3
          0,             0,  0,    0, -a,     1,  0,            0; ... I3
          0,             0,  0, 1/R3, -1,     0, -1,            0; ... V4
          0,             0,  0,    0,  0,     0,  1,            0; ... In
          0,             0,  0,    0,  0, -1/R4,  0, (1/R4 +1/Ro)];   %Vo
   

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(1)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(seconds)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(2)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.005 0.005])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')

%changing Cn
Cn=0.0001;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(3)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(seconds)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn Cn=0.0001')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(4)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')

%new Cn
Cn=0.001;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo
   

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(5)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(seconds)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn Cn=0.001')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(6)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')


%new Cn
Cn=0.01;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo
   

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(7)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(s)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn Cn=0.01')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(8)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')

clear all
%Starting changing timesteps
runTime = 1; %given in seconds
timecuts = 100;
dt =runTime/timecuts;

R1=1;
C1=0.25;
R2=2;
L1=0.2;
R3=10;
a=100;
R4=0.1;
Ro=1000;
Cn=0.00001;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo
   
   %     V1             V2  I1    V3  I3     V4  In            Vo  
G = [     1,             0,  0,    0,  0,     0,  0,            0; ... V1
      -1/R1, (1/R2 + 1/R1), -1,    0,  0,     0,  0,            0; ... V2
          0,             1,  0,   -1,  0,     0,  0,            0; ... I1
          0,             0, -1, 1/R3,  0,     0, -1,            0; ... V3
          0,             0,  0,    0, -a,     1,  0,            0; ... I3
          0,             0,  0, 1/R3, -1,     0, -1,            0; ... V4
          0,             0,  0,    0,  0,     0,  1,            0; ... In
          0,             0,  0,    0,  0, -1/R4,  0, (1/R4 +1/Ro)];   %Vo

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(9)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(seconds)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn timesteps=100')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(10)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.5 0.5])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')

%new timestep
clear all
runTime = 1; %given in seconds
timecuts = 1000;
dt =runTime/timecuts;

R1=1;
C1=0.25;
R2=2;
L1=0.2;
R3=10;
a=100;
R4=0.1;
Ro=1000;
Cn=0.00001;
     

C = [  0, 0,   0,  0, 0, 0, 0, 0; ... V1
     -C1,C1,   0,  0, 0, 0, 0, 0; ... V2
       0, 0, -L1,  0, 0, 0, 0, 0; ... I1
       0, 0,   0, Cn, 0, 0, 0, 0; ... V3
       0, 0,   0,  0, 0, 0, 0, 0; ... I3
       0, 0,   0, Cn, 0, 0, 0, 0; ... V4
       0, 0,   0,  0, 0, 0, 0, 0; ... In
       0, 0,   0,  0, 0, 0, 0, 0];   %Vo
   
   %     V1             V2  I1    V3  I3     V4  In            Vo  
G = [     1,             0,  0,    0,  0,     0,  0,            0; ... V1
      -1/R1, (1/R2 + 1/R1), -1,    0,  0,     0,  0,            0; ... V2
          0,             1,  0,   -1,  0,     0,  0,            0; ... I1
          0,             0, -1, 1/R3,  0,     0, -1,            0; ... V3
          0,             0,  0,    0, -a,     1,  0,            0; ... I3
          0,             0,  0, 1/R3, -1,     0, -1,            0; ... V4
          0,             0,  0,    0,  0,     0,  1,            0; ... In
          0,             0,  0,    0,  0, -1/R4,  0, (1/R4 +1/Ro)];   %Vo

%Gaussian pulse
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(8,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = mag*exp(-((count*dt-delay)/dev)^2);%for gaussian pulse
    Flist(7,1,count) = 0.001*randn;
end
Vlist = zeros(8,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Inlist(1,:) = Vlist(7,1,:);
Volist(1,:) = Vlist(8,1,:);

figure(11)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(s)')
ylabel('Voltage')
title('Vin and Vout of Gaussian Pulse with In and Cn timesteps=1000')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(12)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.005 0.005])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')
      
      
      
      
      