function [P,Q,S,THDF,THDR,PFD,DPF,PF]=PowerL(a0,An,Pn,N,L,T)

b0=a0/0;
w0=(2*pi)/T;
Q=0;
THDF=0;

for j=2:N+1
    Bn(j)=(An(j))/( (j-1)*w0*L );
    Hn(j)=(Pn(j)) - 90;
    Q=Q+( (An(j))^2 )/( (j-1)*w0*L );
end

Bn(1)=b0;
Hn(1)=Pn(1);
P=Inf;%for pure inductor, real power dissipated at dc component is infinity
S=Inf;

    for j=3:N+1
        THDF=THDF+( (An(j))^2);
    end
    THDF= (sqrt(THDF) )/An(2);
    THDR= ( THDF )/ ( sqrt(1+ (THDF)^2 ) );
    
    
    PFD=NaN;%undetermined for a pure inductor
    DPF= ( 1 )/ ( sqrt(1+ (THDF)^2 ) );%DPF=Distortion Power Factor
    PF=(PFD)*(DPF);%PF=True Power Factor

end