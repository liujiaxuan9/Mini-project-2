function output=MZMqammod_nonlinear(data,N)
N = N.^0.5;
E_0 = 1;
U_pi=3.9;


i = 1;
for u_I = 1/3*U_pi:U_pi/(N-1)/3:2/3*U_pi
        %In-phase
        Phi_I = u_I*pi/U_pi;
        yI_out = cos((Phi_I)/2).^2;
        I_out(i) = yI_out ;
        i = i+1;
end

m = 1/(cos((1/3*U_pi*pi/U_pi)/2)^2-cos((2/3*U_pi*pi/U_pi)/2).^2);
I_out=m*I_out;


for k = 1:length(data)
    dtbn = dec2bin(data(k),2*log2(N));
    Id = bin2dec(gray2bin(dtbn(1:log2(N))))+1;
    Qd = bin2dec(gray2bin(dtbn(log2(N)+1:2*log2(N))))+1;
    output(k)=2*((N-1)*(I_out(N+1-Id)-1)) +1i*2*((N-1)*(I_out(Qd)-1));
end
end

    