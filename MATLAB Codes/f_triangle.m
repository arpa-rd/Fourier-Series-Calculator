function [a0,an,bn,An,Pn,F,Wn]=f_triangle(fs,T,amp,xs,vs,N)

w0=(2*pi)/T;
k=-((5*T/2)-1/fs):1/fs:((5*T/2)-1/fs);

syms t   %t is systematic variable now for integrstion
%N=10;
m1=(4*amp)/T;
m2=-m1;

n=1:N;
f1=@(x)   (m1*(x-xs-T/4))+vs  ;
f2=@(x)   (m2*(x-xs-((3*T)/4)))+vs  ;

a0=(1/T)*((int(f1(t),t,xs,xs+T/2))+(int(f2(t),t,xs+T/2,xs+T)));
an=(2/T)*((int(f1(t)*cos(n*w0*(t-xs)),t,xs,xs+T/2))+(int(f2(t)*cos(n*w0*(t-xs)),t,xs+T/2,xs+T)));
bn=(2/T)*((int(f1(t)*sin(n*w0*(t-xs)),t,xs,xs+T/2))+(int(f2(t)*sin(n*w0*(t-xs)),t,xs+T/2,xs+T)));

F=a0;
for i=1:N
   F=F+(an(i)*cos(i*w0*(k-xs)))+(bn(i)*sin(i*w0*(k-xs))); 
end

f=1/T;

Wn=zeros(1,N+1);
An=zeros(1,N+1);
Pn=zeros(1,N+1);
Wn(1)=0;
An(1)=a0;
Pn(1)=0;
for i=2:N+1
Wn(i)=(2*pi*(i-1)*f);
An(i)=sqrt(an(i-1)^2+bn(i-1)^2);
  if an(i-1)~=0
     Pn(i)= -((atan(bn(i-1)/an(i-1)))*(180/pi));
  elseif (an(i-1)==0 && bn(i-1)==0)
     Pn(i)= 0; 
  else
     Pn(i)= -90;
  end
end

end