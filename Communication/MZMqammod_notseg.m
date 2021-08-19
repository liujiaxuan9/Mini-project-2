function output=MZMqammod_notseg(data,N)
N = N.^0.5;
E_0 = 1;
for i = 0:N-1
    I(i+1)= i/(N-1);
end
I = sort(I,'descend');
Phi_ideal = acos(2.*I-1);

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



for i=1:N
    syms x
    x = double(solve(p1*x+p2*x^2 == Phi_ideal(i),x));
    if x(1)>=0
        u(i)=x(1);
    else
        u(i)=x(2);
    end
end



for i=1:N
    %alpha(i)=(a1*(abs(u(i)+U_dc))).^a2;
    alpha(i)=(a1*(abs(u(i)*0.487+U_dc))).^a2;
    I_out(i)= 1/2 * abs(E_0)^2*(1+cos(Phi_ideal(i)))*exp(-alpha(i)*L_total);
end

for k = 1:length(data)
    dtbn = dec2bin(data(k),2*log2(N));
    Id = bin2dec(gray2bin(dtbn(1:log2(N))))+1;
    Qd = bin2dec(gray2bin(dtbn(log2(N)+1:2*log2(N))))+1;
    output(k)=round(2*(N-1)*I_out(N+1-Id)-N+1) +1i*round(2*(N-1)*I_out(Qd)-N+1);
end
end

    