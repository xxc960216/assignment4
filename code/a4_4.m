clear all

disp('model the equation V4 - aI3 - BI3^2 -CI3^3 as a polynomial and solve for it')
disp('This would require a solve for the polynomial per iteration as V4 changes')

runTime = 1; %given in seconds
timecuts = 1000;
dt =runTime/timecuts;

R1=1;
C1=0.25;
R2=2;
L1=0.2;
R3=10;
a=100;
b=50;
c=1;
R4=0.1;
Ro=1000;     

C = [  0, 0,   0, 0, 0, 0, 0; ...
     -C1,C1,   0, 0, 0, 0, 0; ...
       0, 0, -L1, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0];
   
   
   
I3poly = [c b a 0];
I3roots = roots(I3poly);

   
G = [     1,             0,  0,    0,                  0,     0,            0; ...
      -1/R1, (1/R2 + 1/R1), -1,    0,                  0,     0,            0; ...
          0,             1,  0,   -1,                  0,     0,            0; ...
          0,             0, -1, 1/R3,                  0,     0,            0; ...
          0,             0,  0,    0,         I3roots(2),     1,            0; ...
          0,             0,  0, 1/R3,                 -1,     0,            0; ...
          0,             0,  0,    0,                  0, -1/R4, (1/R4 +1/Ro)];
   
      
%Time Stepping function
V1 = 0;      
F = zeros(7,1);
%Flist = [V1; 0; 0; 0; 0; 0; 0];
Flist = zeros(7,1,timecuts);
Flist(1,1,30:timecuts) = 1;
Vlist = zeros(7,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
    Flist(:,:,count) = Vlist(:,:,count);
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Volist(1,:) = Vlist(7,1,:);


figure(1)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(s)')
ylabel('Voltage')
title('Vin and Vout of Step voltage')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(2)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('frequency')
ylabel('magnitude')
title('Fourier transform of output')