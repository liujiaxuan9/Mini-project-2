clear all;
clc;

num = 1000; %Number of bits
M = 256; %Modulation Order 16
data = randi([0 M-1],1,num);
binData = de2bi(data);
modData = MZMqammod(data,M);
refData= qammod(data,M);
o1=real(modData);
o2=imag(modData);

%%%% Oscillator and Balanced Modulator %%%%
fc=1.94*10.^14;              % Carrier Frequency
osr = 128;            %Oversampling rate
fs=osr*fc;            % Sampling Frequency
t=1:10*osr;             % Duration of Signal
s=[];
clear i;
for i=1:1:num       % Modulation (multiplication with carrier signals cos and sin)
    Ac=o1(i);
    As=o2(i);
    s1=Ac*cos(2*pi*(fc/fs)*t);
    s2=As*sin(2*pi*(fc/fs)*t);
    s=[s (s1+s2)];
end
figure(1)
plot(s);
s = awgn(s,10,'measured');
figure(2)
periodogram(s,[],length(s),fs);
%%%% COHERENT DEMODULATION %%%%
ss1=[];
ss2=[];
rs1= 2*cos(2*pi*(fc/fs)*t);  % Sin and Cos generation by Local Oscillator
rs2= 2*sin(2*pi*(fc/fs)*t);
clear i;
for i=1:length(t):length(s)           % Demodulation of Received Signal
ss1=[ss1 rs1.*s(1,i:i+length(t)-1)]; % for I-channel
ss2=[ss2 rs2.*s(1,i:i+length(t)-1)]; % for Q-channel
end
wn=2*(fc/fs);                 
[b1,a1]=butter(2,wn,'low');  

fo1=filtfilt(b1,a1,ss1);       
fo2=filtfilt(b1,a1,ss2);
ro1=[];
ro2=[];
% ro1=resample(fo1,length(o1),length(fo1));
% ro2=resample(fo2,length(o2),length(fo2));
clear i;
for i=1:length(t):length(fo1)           
        ff1=fo1(1,i:i+length(t)-1);
        l1=length(ff1);
        sum1=sum(ff1);
        av1=sum1/l1;
        ro1=[ro1 av1];        
end
clear i;
for i=1:length(t):length(fo2)           
        ff2=fo2(1,i:i+length(t)-1);
        l2=length(ff2);
        sum2=sum(ff2);
        av2=sum2/l2;
        ro2=[ro2 av2];        
end


%%%% Amplitude Detection and Generation of Demodulated Data %%%%
o=o1+1i*o2;
ro=(ro1+1i*ro2);
ro1 = real(ro);
ro2 = imag(ro);
oo1=round(ro1);
oo2=round(ro2);

oo=oo1+1i*oo2;
op=qamdemod(oo,M);
biop = de2bi(op);
figure(3);
subplot(2,1,1);
plot(data,'linewidth',2);
title('Input Data Stream','FontSize',14);

subplot(2,1,2);
plot(op,'linewidth',2);
title('Demodulated Output Data (Same As Input)','FontSize',14);


figure(4)
scatplot = scatterplot((refData),1,0,'r*');
hold on;
scatterplot(ro,1,0,'bx',scatplot);
title('Not segmented MZM');
legend('Original','Received');

[num,ber]=biterr(data,op);
