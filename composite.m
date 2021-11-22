%SRS PLATE
clear all
clc
%Parameters
ax=1; Tinf=0;
lam=1; lam1=0.01;
a1=1.1234; d1=0.7985;
a2=0.0508; d2=0.0492;
a3=0.5088; d3=0.4934;
la1=1; ld1=1; la2=1;
alpha=1.6541;
ba2=0.009; bd2=0.01213;
ta=3.5e6; td=ta;
m=0.02; c11=92776.9;
c12=1.4325; c21=2.0137;
r=1; s=0.01; n=1;
ba1=1; bd1=ba1;
%Coefficients
a14=0; a13=-1/(n*n);
a15d=32*a1*(d1*(r^4)+2*(d2+2*d3)*r*r+1);
a15=-s*s*(a1*(r^4)+1)*(a1-a2*a2)*c11/a15d;
a16d=pi*pi*(d1*(r^4)+2*(d2+2*d3)*r*r+1);
a16=(bd1*r*r+1)*c12/a16d;
a17d=3*(pi^4)*(a3+r*r*(-a2*(a2+2*a3)+a1*(1+r*r*a3)))*(1+(r^4)*d1+2*r*r*(d2+2*d3));
a17=32*r*r*s*s*a3*c11*c21*(-a2+ba1+r*r*(a1-a2*ba1))/a17d;
a24d=s*s*n*((a3+r*r*(-a2*(a2+2*a3)+a1*(1+r*r*a3)))*ta+c21*(a3+r*r*(a1+ba1*(-2*a2+ba1+a3*(-2+r*r*ba1))))*ba2*ta);
a24=(a3+r*r*(-a2*(a2+2*a3)+a1*(1+r*r*a3)))*(pi*pi*s*s*(1+r*r*la1)-4*m*la2)/a24d;
a25dd=a3+r*r*(-a2*(a2+2*a3)+a1*(1+r*r*a3));
a25d=3*a1*((c21*ba2*ta*(a3+r*r*(a1+ba1*(-2*a2+ba1+a3*(-2+r*r*ba1))))/a25dd)+ta);
a25=-ba2*4*(a1*(-1+r*r*(a2-ba1))+a2*ba1)*ta/a25d;
a26dd=(a1*(a3*r*r+1)-a2*(a2+2*a3))*r*r+a3;
a26d=pi*pi*(s^4)*((c21*((a1+ba1*(-2*a2+ba1+a3*(r*r*ba1-2)))*r*r+a3)/a26dd)+ta);
a26=64*m*la2*alpha/a26d;
a32=-(pi*pi*s*s*(r*r*la1+1)+12*(m+1)*la2)/(n*s*s*ta);
a33=(1)*bd2*(1+r*r*bd1)*pi*pi;

%Composite to ortho
Ak1=-(a13+a14); Ak2=-a16; Ak3=-a15; Ak4=a17;
Bk1=a32; Bk2=-a33;
Ck1=a24; Ck2=-a26; Ck3=-a25;

%Initial Condition
A=[1]; B=[0]; C=[0]; y=[0];
h=pi/300; N=90000/pi;
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
plot(time,A)
hold on
%pks=findpeaks(A)
b31=-a33*a16/a32