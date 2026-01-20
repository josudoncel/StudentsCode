function [Vberria,politika]= Value_iteration(E,P,R,lambda)

% ABIARAPENAK
V = ones(length(E),1);
Vberria = zeros(length(E),1);
politika = zeros(length(E),1);

% BUKLE NAGUSIA
while norm(V-Vberria,inf)>10^(-10)
    V = Vberria;
    for e = 1:length(E)
        V1 = R(e,1)+lambda*(P(e,1)*V(1,1)+(1-P(e,1))*V(2,1)); %1. akzioa (Txokolaterik ez jan) hartzean itxarotako balioa
        V2 = R(e,2)+lambda*(P(e,2)*V(1,1)+(1-P(e,2))*V(2,1)); %2. akzioa (Txokolatea jan) hartzean itxarotako balioa
        Vberria(e)=max(V1,V2); %egoera bakoitzeko balio optimoak
        if V1 > V2
            politika(e,1)=1; %1. akzioa hobea (Txokolaterik ez jan)
        elseif
            politika(e,1)=2; %2. akzioa hobea (Txokolatea jan)
        else
            politika(e,1)=3; %3. akzioa hobea (Berdin dio)
        end
    end
end
