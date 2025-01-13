function [a0,an,bn,An,Pn,F,Wn]=f_rect(fs,T,W,amp,xs,vs,N)

w0=(2*pi)/T;
k=-(2*T+(W/2)-xs)-1/fs:1/fs:(2*T+(W/2)+xs)+1/fs;

syms t   %t is systematic variable now for integrstion
n=1:N;

a0=(1/T)*((int(amp+vs,t,-(W/2-xs),W/2+xs))+(int(0+vs,t,(W/2+xs),(T-W/2+xs))));
an=(2/T)*((int((amp+vs)*cos(n*w0*t),t,-(W/2-xs),W/2+xs))+(int((0+vs)*cos(n*w0*t),t,(W/2+xs),(T-W/2+xs))));
bn=(2/T)*((int((amp+vs)*sin(n*w0*t),t,-(W/2-xs),W/2+xs))+(int((0+vs)*sin(n*w0*t),t,(W/2+xs),(T-W/2+xs))));

F=a0;
for i=1:N
   F=F+(an(i)*cos(i*w0*k))+(bn(i)*sin(i*w0*k)); 
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