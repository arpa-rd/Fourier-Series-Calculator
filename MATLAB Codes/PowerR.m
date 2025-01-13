function [P,Q,S,THDF,THDR,PFD,DPF,PF]=PowerR(a0,An,Pn,N,R)

b0=a0/R;
P=((a0)^2)/R;
THDF=0;

for j=2:N+1% Here An(1)=a0 don't be confused
    Bn(j)=An(j)/R;%Bn,Hn has N+1 term starting from 1 to N+1
    Hn(j)=Pn(j);
    P=P+( (An(j))^2 )/R;
end
Bn(1)=b0;
Hn(1)=Pn(1);

    Q=0;
    S=P;
    
    for j=3:N+1
        THDF=THDF+( (An(j))^2);
    end
    THDF= (sqrt(THDF) )/An(2);
    THDR= ( THDF )/ ( sqrt(1+ (THDF)^2 ) );
    
    PFD=1;%PFD=Displacement Power Factor
    DPF= ( 1 )/ ( sqrt(1+ (THDF)^2 ) );%DPF=Distortion Power Factor
    PF=(PFD)*(DPF);%PF=True Power Factor

end