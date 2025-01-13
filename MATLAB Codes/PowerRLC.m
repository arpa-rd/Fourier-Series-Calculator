function [P,Q,S,THDF,THDR,PFD,DPF,PF]=PowerRLC(a0,An,Pn,N,R,L,C,T)

b0=a0/Inf;
w0=(2*pi)/T;
P=0;
Q=0;
THDF=0;

for j=2:N+1
    Zeq=   R+( (((j-1)*w0*L)- (1/( (j-1)*w0*C) ) )*i  );% real and imag
    Z=abs(Zeq);
    ang=(angle(Zeq))*(180/pi);
    
    Bn(j)= An(j)/Z;
    Hn(j)= Pn(j)-ang;
end
Bn(1)=b0;
Hn(1)=Pn(1);

for j=1:N+1
    P=P+  An(j)*Bn(j)*cos(  ((Pn(j)-Hn(j))*(pi/180))   ) ;
    Q=Q+ An(j)*Bn(j)*sin (  ((Pn(j)-Hn(j))*(pi/180))   ) ;
end
S=sqrt( (P^2)+ (Q^2) );


    for j=3:N+1
        THDF=THDF+( (An(j))^2);
    end
    THDF= (sqrt(THDF) )/An(2);
    THDR= ( THDF )/ ( sqrt(1+ (THDF)^2 ) );
    
    PFD=(P/S);%PFD=Displacement Power Factor
    DPF= ( 1 )/ ( sqrt(1+ (THDF)^2 ) );%DPF=Distortion Power Factor
    PF=(PFD)*(DPF);%PF=True Power Factor

end