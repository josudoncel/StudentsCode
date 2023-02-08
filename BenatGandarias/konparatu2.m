function konparatu2

l=0.6;
B=5;
mu=1; 
r_FIFO=FIFO(B,l,mu);
r_PS=ProcessorSharing(B,l,mu);
end


function [Pi]=Banaketak(B,a,m)
G=zeros(B+1);
G(1,1)=-a; G(B+1,B+1)=-m; G(1,2)=m; G(B+1,B)=a;
for i=2:B
    G(i,i-1)=a;
    G(i,i+1)=m;
    G(i,i)=-(a+m);
end

G(B+2,1:B+1)=1;
S=zeros(B+1,1);
S(B+2,1)=1;
Pi=linsolve(G,S);
end

function [Itxaropena_AoI]=FIFO(B,a,m)
C=sym(zeros(2*B+1,B+1));
V=sym('v',[B+1,B+1]);

%C matrizea sortzen
C(1,1)=V(1,1); C(2,1)=V(2,2);
%Errenkada bikoitiak sortzen
batura=1;
for i=4:2:2*B
   C(i,1:i/2)=V(batura+2,2:batura+2);
   batura=batura+1;
end
%Errenkada bakoitiak sortzen
suma=2;
for i=3:2:2*B+1
    C(i,1:suma)=V(suma,1:suma);
    suma=suma+1;
end

%Ekuazio sistemako eskuineko aldeko matrizea sortzen
K=sym(zeros(B+1,B+1));
K(1,1:B+1)=m*C(2,1:B+1);
batura2=1;
for i=1:B-1
    K(i+1,1:B+1)=a*C(batura2,1:B+1)+m*C(batura2+3,1:B+1);
    batura2=batura2+2;
end
K(B+1,1:B+1)=a*C(2*B-1,1:B+1)+a*C(2*B+1,1:B+1);

[Pi]=Banaketak(B,a,m);
Pisuak=zeros(B+1,B+1);
for i=1:B+1
    for j=1:i
        Pisuak(i,j)=1;
    end
end
EskuinekoMat=Pisuak.*Pi+K;
%Ekuazio sistemako ezkerreko aldeko matrizea sortzen
t=zeros(B+1,1);
t(1)=a;
for i=2:B+1
    t(i)=a+m;
end
EzkerrekoMat=V.*t;

%Ekuazio sistema ebazten
kontatzaile=1;
for i=1:B+1
    for j=1:B+1
    eqns(kontatzaile)=[EzkerrekoMat(i,j)==EskuinekoMat(i,j)];
    kontatzaile=kontatzaile+1;
    end
end
Aldagaiak=solve(eqns,V);

%Age of Informationaren itxaropen matematikoa kalkulatzen
fns=fieldnames(Aldagaiak);
Itxaropena_AoI=0;
for i=1:B+1
    Itxaropena_AoI=Itxaropena_AoI+Aldagaiak.(fns{i});
end
Itxaropena_AoI=vpa(Itxaropena_AoI);

end

function [Itxaropena_AoI]=ProcessorSharing(B,a,m)
n=0;
if B==1
    n=0;
end
if B>=2
    for i=0:B-2
        n=n+3+i;
    end
end
C=sym(zeros(3+n,B+1));
V=sym('v',[B+1,B+1]);

%C matrizea sortzen
C(1,1)=V(1,1); C(2,1)=V(2,2);C(3,1:2)=V(2,1:2);
%Lambda tasa duten errenkadak sortzen
errenkada=0;
if B>=2
    for i=0:B-2
        errenkada=errenkada+3+i;
        C(3+errenkada,1:i+3)=V(i+3,1:i+3);
    end
end
%Probabilitate banatua(m/n modukoa) duten errenkadak sortzen    
if B>=2
    errenk=4;
    d=0;
    for i=0:B-2
        errenk=errenk+d;
        d=3+i;
        for j=0:i+1
            C(errenk+j,1:i+2)=V(i+3,2:i+3);
        end
    end
    %Gure aldagai imajinarioak sartzen
    errenka=5;
    de=0;
    for i=0:B-2
        errenka=errenka+de;
        de=3+i;
        for j=0:i
            C(errenka+j,1:j+1)=C(errenka+j,j+2);
        end
    end
end


%Ekuazio sistemako eskuineko aldeko matrizea sortzen
K=sym(zeros(B+1,B+1));
K(1,1:B+1)=m*C(2,1:B+1);
K(B+1,1:B+1)=a*C(3+n-(B+1),1:B+1)+a*C(3+n,1:B+1);
erren=4;
s=0;
batura2=0;
sum=2;
for i=0:B-2
    erren=erren+s;
    s=3+i;
    batura=sym(zeros(1,B+1));
    for j=0:i+1
        batura=batura+(m/(i+2))*C(erren+j,1:B+1);
    end
    
    
    batura2=-1+batura2+sum;
    
    K(i+2,1:B+1)=a*C(batura2,1:B+1)+batura;
    sum=sum+1;
end


[Pi]=Banaketak(B,a,m);
Pisuak=zeros(B+1,B+1);
for i=1:B+1
    for j=1:i
        Pisuak(i,j)=1;
    end
end
EskuinekoMat=Pisuak.*Pi+K;
%Ekuazio sistemako ezkerreko aldeko matrizea sortzen
t=zeros(B+1,1);
t(1)=a;
for i=2:B+1
    t(i)=a+m;
end
EzkerrekoMat=V.*t;

%Ekuazio sistema ebazten
kontatzaile=1;
for i=1:B+1
    for j=1:B+1
    eqns(kontatzaile)=[EzkerrekoMat(i,j)==EskuinekoMat(i,j)];
    kontatzaile=kontatzaile+1;
    end
end
Aldagaiak=solve(eqns,V);

%Age of Informationaren itxaropen matematikoa kalkulatzen
fns=fieldnames(Aldagaiak);
Itxaropena_AoI=0;
for i=1:B+1
    Itxaropena_AoI=Itxaropena_AoI+Aldagaiak.(fns{i});
end
Itxaropena_AoI=vpa(Itxaropena_AoI);

end
