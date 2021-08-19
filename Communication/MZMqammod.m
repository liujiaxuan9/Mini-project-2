function output=MZMqammod(data,N)
N = N.^0.5;
E_0 = 1;
% for i = 0:N-1
%     I(i+1)= i/(N-1);
% end
% I = sort(I,'descend');
I = linspace(0,1,N);
Phi_ideal = acos(2.*I-1);
% Phi_difference(1)=Phi_ideal(1);
% for i=2:N
%     Phi_difference(i)= Phi_ideal(i)-Phi_ideal(i-1);
% end
Phi_difference= Phi_ideal(2:end)-Phi_ideal(1:end-1);

U_pi = 3.9;
%p1=1.8;
%p2=0.03305;

p1= -0.2856*U_pi+ 1.6051;
p2 = 0.0817;

U_dc = -10;
L_total=0.0006;
%a1=0.2079;
%a2=4.5285;
a1 = -0.0699*U_pi+ 0.47;
a2 = 2.395*U_pi - 4.4528;

for i = 1:N/2+1
    Phi(i)= Phi_ideal(i);
    %u(i)= Phi_difference(i)/pi*U_pi;
    %u(i)= solve(p1*x+p2*x^2 == Phi(i),x);
    %I_out(i)= 1/2 * abs(E_0)^2*(1+cos(Phi(i)));
end
%u(N/2+1)=Phi(N/2+1)/pi*U_pi;
%u(N/2+1)=solve('p1*u(N/2+1)+p2*u(N/2+1)^2=Phi(N/2+1)','u(i)');


for i=1:N/2+1
    syms x
    x = double(solve(p1*x+p2*x^2 == Phi(i),x));
    if x(1)>=0
        u(i)=x(1);
    else
        u(i)=x(2);
    end
end

u1(1)=u(1);
%Phi1(1)=Phi(1);
%Phi1_difference(1)=Phi1(1);
for i=2:N/2
    u1(i)=u(i)-u(i-1);
    %Phi1(i) = Phi(i-1)+p1*u1(i)+p2*u1(i)^2;
    %Phi1_difference(i)=Phi1(i)-Phi1(i-1);
end
u1(N/2+1)=u(N/2+1);
%Phi1(N/2+1)=p1*u1(N/2+1)+p2*u1(N/2+1)^2;
%Phi1_difference(N/2+1)=Phi1(N/2+1);


for i=1:N/2+1
    %alpha(i)=(a1*(abs(u(i)+U_dc))).^a2;
    alpha(i)=(a1*(abs(u(i)*0.487+U_dc))).^a2;
    L(i)=u(i)/U_pi*L_total;
    I_out(i)= 1/2 * abs(E_0)^2*(1+cos(Phi(i)))*exp(-alpha(i)*L(i));
end


    
a=1;
for i = N/2+2:N
    Phi(i)= Phi(i-1)+Phi_difference(N/2+1-a);
    u1(i)=u1(i-1)+u1(N/2+1-a);
    %alpha(i)=(a1*(abs(u1(i)+U_dc))).^a2;
    alpha(i)=(a1*(abs(u1(i)*0.487+U_dc))).^a2;
    L(i)=u1(i)/U_pi*L_total;
    %L(i)=u1(i)/5.4852*L_total;
    %Phi1(i)=Phi1(i-1)+Phi1_difference(N/2+1-a);
    I_out(i)= 1/2 * abs(E_0)^2*(1+cos(Phi(i)))*exp(-alpha(i)*L(i));
    a = a+1;
end

for k = 1:length(data)
    dtbn = dec2bin(data(k),2*log2(N));
    Id = bin2dec(gray2bin(dtbn(1:log2(N))))+1;
    Qd = bin2dec(gray2bin(dtbn(log2(N)+1:2*log2(N))))+1;
    output(k)=round(2*(N-1)*I_out(N+1-Id)-N+1) +1i*round(2*(N-1)*I_out(Qd)-N+1);
end
end

    