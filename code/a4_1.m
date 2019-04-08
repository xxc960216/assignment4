clear all


R1=1;
C1=0.25;
R2=2;
L1=0.2;
R3=10;
a=100;
R4=0.1;
Ro=1000;


%       V1    V2     V3          V5      IL3
% G = [-1/R1, 0,     0,          0,      0; ...%N1
%      1/R1 , -1/R2, 0,          0,      0; ...%N2
%      0,     0,     -1/R3,      0,      0; ...%N3
%      0,     0,     -a/(R3*R4), -1/R4,  0; ...%N4
     

C = [  0, 0,   0, 0, 0, 0, 0; ...
     -C1,C1,   0, 0, 0, 0, 0; ...
       0, 0, -L1, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0];
   
G = [     1,             0,  0,    0,  0,     0,            0; ...
      -1/R1, (1/R2 + 1/R1), -1,    0,  0,     0,            0; ...
          0,             1,  0,   -1,  0,     0,            0; ...
          0,             0, -1, 1/R3,  0,     0,            0; ...
          0,             0,  0,    0, -a,     1,            0; ...
          0,             0,  0, 1/R3, -1,     0,            0; ...
          0,             0,  0,    0,  0, -1/R4, (1/R4 +1/Ro)];

      
V1 = 10;      
F = [V1; 0; 0; 0; 0; 0; 0];

w = 0;
V = (G+1i*w*C)\F;

for k =1:21
    vp = -10 +k -1;
    F(1,1) = vp;
    
    V(:,:,k) = (G+1i*w*C)\F;
end

Vo(1,:) = V(7,1,:);
V3(1,:) = V(4,1,:).*R3;

figure(1)
plot(-10:1:10,Vo)
title('DC: Vo for -10 to 10 V V1')

figure(2)
plot(-10:1:10,V3)
title('DC: V3 for -10 to 10 V V1')


%changing omega
F(1,1) =10;
for w = 1:1000
    V(:,:,w) = (G+1i*w*C)\F;
end

clear Vo
Vo(1,:) = V(7,1,:);

Vo1 = 20*log10(Vo/V1);

figure(3)
semilogx(1:1000,Vo1)
title('AC: Gain in dB with varying w')

%Random perturbations on C's
w=pi;
std = 0.05;
for i = 1:100
    Cnew = normrnd(C1,std);
    C(2,1) = -Cnew;
    C(2,2) = Cnew;
    
    V(:,:,i) = (G+1i*w*C)\F;
end

clear Vo
Vo(1,:) = V(7,1,:);

Vo1 = 20*log10(Vo/V1);

figure(4)
hist(real(Vo1(:)))
title('AC: Gain with varying C')


