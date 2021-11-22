clear all
clc

%Material Properties
pk12=1;
v12=0.064; v21=pk12*v12;
ax=0.46*6.6739; %alphax=3.07alpha
ay=3.07*6.6739; %alphay=0.46*alpha
mu=1-v21*v12;
px=0; py=0; %Body force

%Parameter Values
n=1; c1=0.164; c2=0.256;
lam=1; lam1=0.01;
k1=1;  k2=1;
f=0;   Bb=0.13; B1=1;%Byy/Bxx;
tC=3.5e6;
M=0.02;
Tinf=0;

%Constants
s2=1;%(b*b)/(a*a); 
s4=1;%(b^4)/(a^4);
%a1=(pi^4)*b/(4*(a^3)); a2=(pi^4)*a/(4*(b^3)); a3=(pi^4)/(4*a*b);
a4=-(pi^2)/(4*lam);    a5=-lam*(pi^2)/(4);    %a7=(4*a*b)/(pi^2);
%a8=(a*b)/4;            c1=(mu*G12)/E1;       c2=(E2/E1);
c3=v21+2*c1;           l1=ax+v21*ay; 
l2=v21*ax+c2*ay; %Dx=(E1*(h^3))/(12*mu);

X1=a4+k2*a5+M/(lam1*lam*lam);%a4+(kyy*a5/kxx)+(2*H*a8/(kxx*h));
X2=(0.25/lam)*(1+Bb+(ay*B1*Bb/ax));
X3=3*(1+(ay/ax))/(1+(1/c2)+((mu/c1)-2*v12));%(3*s2*E2*G12*(1+s2*ay/ax)/(s4*E1*G12+E2*G12+s2*E2*(E1-2*G12*v12)));
X4=3*(1+(ay/ax))/(1+c2+c2*((mu/c1)-2*v12));%(3*s2*E1*G12*(1+s2*ay/ax)/(s4*E1*G12+E2*G12+s2*E2*(E1-2*G12*v12)));
%g2=(Dx*(pi^4))/(rho*h*a^4)*(1+c2*(a^4/b^4)+2*c3*(a*a/b*b))-(pi^2/rho*a*b*h)*((px*b/a)+(py*a/b));
%g3=(4*Dx*a/(rho*(a^4)*b*h))*((l1*a4/ax)+(l2*a5/ax));
%g4=((3*mu*Dx*(pi^4)*a*a)/(8*rho*h*(a^4)*b*b))*((a*a*E2)/(b*b*E1)+b*b/a*a);
%g5=((Dx*128*mu)/(rho*h*(a^4)))*(1+s2*ay/ax)/(s4*E1/E2+1+s2*((E1/G12)-2*v12));
%g7=(kxx/(C*a8))*(a4+a5*kyy/kxx+(12*kzz*a8*(M-1))/(kxx*h*h));
%g8=((a*a*ax*T0*Bxx)/(a8*C))*(a4+Byy*a5/Bxx);
%g9=(-a*a*C*X2/kxx)+(a*a*ax*T0*Bxx/(12*kxx))*((1-(v12*Byy/Bxx))*(a*X3/b)+((Byy/Bxx)-v21)*b*X4/a);
%g10=X1; 
%g11=(kzz*M*a*a*4*a7)/(kxx*(h^4));
%g12=((a*a*C*ax*T0*Bxx)/(12*kxx*C))*((1-(v12*Byy/Bxx))*(2*b/a)+((Byy/Bxx)-v21)*2*a/b);

%Non Dimensionalised parameters
pb=0;%(((px*b/a)+(py*a/b))*(a^4))/(Dx*a*b*pi*pi*(1+c2/s4+2*c3/s2));
Ak1=(1-pb)/(n*n);
Ak2=4*lam*((l1*a4/ax)+(l2*a5/ax))/(n*n*(pi^4)*(1+c2*(lam^4)+2*c3*lam*lam));
Ak3=3*mu*(c2*(lam^4)+1)/(8*n*n*(1+c2*(lam^4)+2*c3*lam*lam));
Ak4=128*mu*(1+(ay/lam*lam*ax))/(n*n*(pi^4)*((1/c2)+1+(1/lam^2)*(mu/c1-2*v12))*(1+c2*(lam^4)+3*c3*lam*lam));
%fb=(16*f*(a^4))/((pi^6)*Dx*h*(1+c2*(lam^4)+2*c3*lam*lam));
Bk1=(4*lam/(n*tC))*(a4+a5*k2-12*k1*(M+1)*(0.25/(lam1*lam1*lam)));
Bk2=(Bb*4*lam)*(a4+B1*a5);
X5=-X2+(Bb/12)*((1-v12*B1)*lam*X3+(B1-v21)*X4/lam);
Ck1=(X1/(n*tC*X5));
Ck2=k1*M*16/(n*tC*X5*(lam1^4)*(pi^2)*lam);
Ck3=(Bb/(3*X5))*((1-v12)*(B1/lam)+(B1-v21)*lam);

%Initial Condition
A=[1]; B=[0]; C=[0]; y=[0];
h=pi/300; N=60000/pi;
time=0:h:N*h;
for i=1:N
        k1=e(y(i));                                m1=fn(A(i),B(i),C(i),Ak1,Ak2,Ak3,Ak4);                            
        n1=g(B(i),y(i),Bk1,Bk2);                   o1=hn(A(i),C(i),Ck1,Ck2,Ck3,ax,Tinf,y(i));
        k2=e(y(i)+m1*h/2);                         m2=fn(A(i)+k1*0.5*h,B(i)+n1*0.5*h,C(i)+o1*0.5*h,Ak1,Ak2,Ak3,Ak4); 
        n2=g(B(i)+n1*0.5*h,y(i)+m1*0.5*h,Bk1,Bk2); o2=hn(A(i)+k1*0.5*h,C(i)+o1*0.5*h,Ck1,Ck2,Ck3,ax,Tinf,y(i)+m1*0.5*h);
        k3=e(y(i)+m2*h/2);                         m3=fn(A(i)+k2*0.5*h,B(i)+n2*0.5*h,C(i)+o2*0.5*h,Ak1,Ak2,Ak3,Ak4); 
        n3=g(B(i)+n2*0.5*h,y(i)+m2*0.5*h,Bk1,Bk2); o3=hn(A(i)+k2*0.5*h,C(i)+o2*0.5*h,Ck1,Ck2,Ck3,ax,Tinf,y(i)+m2*0.5*h);
        k4=e(y(i)+m3*h);                           m4=fn(A(i)+k3*h,B(i)+n3*h,C(i)+o3*h,Ak1,Ak2,Ak3,Ak4);             
        n4=g(B(i)+n3*h,y(i)+m3*h,Bk1,Bk2);         o4=hn(A(i)+k3*h,C(i)+o3*h,Ck1,Ck2,Ck3,ax,Tinf,y(i)+m3*h);
        A(i+1)=A(i)+(k1+2*k2+2*k3+k4)*h/6;
        y(i+1)=y(i)+(m1+2*m2+2*m3+m4)*h/6;
        B(i+1)=B(i)+(n1+2*n2+2*n3+n4)*h/6;
        C(i+1)=C(i)+(o1+2*o2+2*o3+o4)*h/6;
end
%A=A(1,1500:end);
%time=time(1,1500:end);
%plot(time,A,'-o')
[pks,locs]=findpeaks(A);
locs=locs*h;
plot(locs,pks,'-o')
hold on