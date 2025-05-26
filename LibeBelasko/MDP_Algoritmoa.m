function [Vberria,politika]= MDP_Algoritmoa(N,cv,ci0,ci1,lambda,gamma,rho,alpha,delta1,delta2)
% Hipotesia: N(gamma+alpha+delta1+delta2)<1
% Egoeren multzoa definitu.
E=[];
for i=0:N
    for j=0:N
        for l=0:N
            if i+j+l<=N
                E=[E;[i,j,l]];
            end
        end 
    end
end
[n,~]=size(E);
% V, V^*, politika eta kontadore bat abiarazi.
V=ones(n,1); 
Vberria=zeros(n,1);
politika = zeros(n,1);
k=0; % Kontadorea probabilitate gustiak ondo definituta daudela egiaztatzeko erabiliko dugu.
% Konbergentziarako buklea
while norm(V-Vberria,inf) > 10^(-4)
        V=Vberria;
        % Egoera bakoitzaren arabera Bellmanen ekuazioa erabiliz egoera
        % balio funtzio optimoa kalkulatuko dugu.
        for i = 1:n
            % (0,0,0) egoera.
            if all(E(i,:)==zeros(3,1))
                [~, errenkada1] = ismember([1,0,0], E, 'rows');
                [~, errenkada2] = ismember([0,0,0], E, 'rows');
                Vberria(i)= lambda*(nchoosek(N,1)*alpha*(1-alpha)^(N-1)*V(errenkada1)+(1-nchoosek(N,1)*alpha*(1-alpha)^(N-1))*V(errenkada2));
                politika(i)=0;
                % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                if nchoosek(N,1)*alpha*(1-alpha)^(N-1) <0 || nchoosek(N,1)*alpha*(1-alpha)^(N-1)>1|| 1-nchoosek(N,1)*alpha*(1-alpha)^(N-1)<0 || 1-nchoosek(N,1)*alpha*(1-alpha)^(N-1)>1
                    k=k+1;
                end
            end
            % (i,0,0) egoera non 0<i<=N.
            if E(i,1)~=0 && E(i,2)== 0 && E(i,3)== 0
                [~, errenkada1] = ismember([E(i,1)-1,0,0], E, 'rows');
                [~, errenkada2] = ismember([E(i,1)+1,0,0], E, 'rows');
                [~, errenkada3] = ismember([E(i,1),0,0], E, 'rows');
                % i=1 bada.
                if E(i,1)==1
                    pi=0;
                    V0= (E(i,1)/N)*cv*pi+lambda*(pi*(1-alpha)^(N-1)*V(errenkada1)+nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)*(1-pi)*V(errenkada2)...
                    +(1-pi*(1-alpha)^(N-1)-nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)*(1-pi))*V(errenkada3));

                    pi=1;
                    V1= (E(i,1)/N)*cv*pi+lambda*(pi*(1-alpha)^(N-1)*V(errenkada1)+nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)*(1-pi)*V(errenkada2)...
                    +(1-pi*(1-alpha)^(N-1)-nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)*(1-pi))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if (1-alpha)^(N-1)<0 || (1-alpha)^(N-1)>1 || nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)<0 || nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)>1 || 1-(1-alpha)^(N-1)<0 || 1-(1-alpha)^(N-1)>1 || 1-nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)<0 || 1-nchoosek(N-1,1)*alpha*(1-alpha)^(N-2)>1
                        k=k+1;
                    end
                % i=N bada.
                elseif E(i,1)==N
                    pi=0;
                    V0= (E(i,1)/N)*cv*pi+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*V(errenkada1)...
                        +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1))*V(errenkada3));

                    pi=1;
                    V1= (E(i,1)/N)*cv*pi+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*V(errenkada1)...
                        +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1))*V(errenkada3));
                % 1<i<N bada.
                else
                    pi=0;
                    V0= (E(i,1)/N)*cv*pi+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-alpha)^(N-E(i,1))*V(errenkada1)...
                        +nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)*(1-pi)^(E(i,1))*V(errenkada2)...
                        +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-alpha)^(N-E(i,1))-nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)*(1-pi)^(E(i,1)))*V(errenkada3));

                    pi=1;
                    V1= (E(i,1)/N)*cv*pi+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-alpha)^(N-E(i,1))*V(errenkada1)...
                        +nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)*(1-pi)^(E(i,1))*V(errenkada2)...
                        +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-alpha)^(N-E(i,1))-nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)*(1-pi)^(E(i,1)))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)<0 || nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)>1 || 1-nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)<0 || 1-nchoosek(N-E(i,1),1)*alpha*(1-alpha)^(N-E(i,1)-1)>1
                        k=k+1;
                    end
                end
                % Minimoa aukeratu.
                if V0<V1
                    Vberria(i)=V0;
                    politika(i)=0;
                else
                    Vberria(i)=V1;
                    politika(i)=1;
                end
            end
            
            % (0,j,0) egoera non 0<j<=N.
            if E(i,1)== 0 && E(i,2)~=0 && E(i,3)== 0
                [~, errenkada1] = ismember([0,E(i,2)-1,0], E, 'rows');
                [~, errenkada2] = ismember([1,E(i,2),0], E, 'rows');
                [~, errenkada3] = ismember([0,E(i,2),0], E, 'rows');
                % j=N bada.
                if E(i,2)==N
                    Vberria(i)=ci0*(E(i,2)/N)+lambda*(nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*V(errenkada1)...
                                +(1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if  nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)<0 ||nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)>1 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)<0 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)>1
                        k=k+1;
                    end 
                % j<N bada.
                else
                    Vberria(i)=ci0*(E(i,2)/N)+lambda*(nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))*V(errenkada1)...
                                +nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2))*V(errenkada2)...
                                +(1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))-nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2)))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if  nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))<0 ||nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))>1||nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2))<0 ||nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2))>1|| 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))-nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2))<0 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-alpha)^(N-E(i,2))-nchoosek(N-E(i,2),1)*alpha*(1-alpha)^(N-E(i,2)-1)*(1-delta1)^(E(i,2))>1
                        k=k+1;
                    end 
                end
                politika(i)=0;
            end
            % (0,0,k) egoera non 0<k<=N.
            if E(i,1)== 0 && E(i,2)==0 && E(i,3)~= 0
                [~, errenkada1] = ismember([0,0,E(i,3)-1], E, 'rows');
                [~, errenkada2] = ismember([1,0,E(i,3)], E, 'rows');
                [~, errenkada3] = ismember([0,0,E(i,3)], E, 'rows');
                %k=N bada.
                if E(i,3)==N
                    Vberria(i)=ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*V(errenkada1)...
                                +(1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if  nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 ||nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1
                        k=k+1;
                    end 
                % k<N bada.
                else
                    Vberria(i)=ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))*V(errenkada1)...
                                +nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3))*V(errenkada2)...
                                +(1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))-nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3)))*V(errenkada3));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if  nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))<0 ||nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))>1||nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3))<0 ||nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3))>1|| 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))-nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,3))-nchoosek(N-E(i,3),1)*alpha*(1-alpha)^(N-E(i,3)-1)*(1-delta2)^(E(i,3))>1
                        k=k+1;
                    end 
                end
                politika(i)=0;
            end
            % (i,j,0) egoera non i+j<=N eta i,j>=1.
            if E(i,1)~=0 && E(i,2)~=0 && E(i,3)== 0
                [~, errenkada1] = ismember([E(i,1)-1,E(i,2),0], E, 'rows');
                [~, errenkada2] = ismember([E(i,1)-1,E(i,2)+1,0], E, 'rows');
                [~, errenkada3] = ismember([E(i,1)-1,E(i,2),1], E, 'rows');
                [~, errenkada4] = ismember([E(i,1),E(i,2)-1,0], E, 'rows');
                [~, errenkada5] = ismember([E(i,1)+1,E(i,2),0], E, 'rows');
                [~, errenkada6] = ismember([E(i,1),E(i,2),0], E, 'rows');
                % i+j=N bada.
                if E(i,1)+E(i,2)==N
                    %i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(pi*(1-delta1)^(E(i,2))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*V(errenkada4)...
                            +(1-pi*(1-delta1)^(E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1))))*V(errenkada6));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(pi*(1-delta1)^(E(i,2))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*V(errenkada4)...
                            +(1-pi*(1-delta1)^(E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1))))*V(errenkada6));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if (1-delta1)^(E(i,2))<0 || (1-delta1)^(E(i,2))>1 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))<0 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))>1 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))<0 ||gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))>1 || 1-(1-delta1)^(E(i,2))<0 || 1-(1-delta1)^(E(i,2))>1 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))<0 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))>1
                            k=k+1;
                        end
                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*V(errenkada4)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)))*V(errenkada6));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*V(errenkada4)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)))*V(errenkada6));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))>1 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))>1 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))<0 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))>1
                            k=k+1;
                        end
                    
                    end
                % i+j<N bada.
                else
                    % i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(pi*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada4)...
                            +nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2)))*V(errenkada6));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(pi*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada4)...
                            +nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2)))*V(errenkada6));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if (1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || (1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))<0 || nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))>1|| 1-(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || 1-(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))<0 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))>1
                            k=k+1;
                        end
                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada4)...
                            +nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2)))*V(errenkada6));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ ci0*(E(i,2)/N)+ lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))*V(errenkada4)...
                            +nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2)))*V(errenkada6));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))>1 || nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))<0 || nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))>1 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))<0 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,2))-nchoosek(N-E(i,1)-E(i,2),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))>1
                            k=k+1;
                        end
                    end
                end
                if V0<V1
                    Vberria(i)=V0;
                    politika(i)=0;
                else
                    Vberria(i)=V1;
                    politika(i)=1;
                end
            end
            
            % (i,0,k) egoera non i+k<=N eta i,k>=1.
            if E(i,1)~=0 && E(i,2)==0 && E(i,3)~=0
                [~, errenkada1] = ismember([E(i,1)-1,0,E(i,3)], E, 'rows');
                [~, errenkada2] = ismember([E(i,1),0,E(i,3)-1], E, 'rows');
                [~, errenkada3] = ismember([E(i,1)+1,0,E(i,3)], E, 'rows');
                [~, errenkada4] = ismember([E(i,1),0,E(i,3)], E, 'rows');
                % i+k=N bada.
                if E(i,1)+E(i,3)==N
                    % i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci1 *(E(i,3)/N)+lambda*(pi*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*V(errenkada2)...
                            +(1-pi*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi))*V(errenkada4));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci1 *(E(i,3)/N)+lambda*(pi*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*V(errenkada2)...
                            +(1-pi*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi))*V(errenkada4));
                        
                        if (1-delta2)^(E(i,3))<0|| (1-delta2)^(E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1 || 1-(1-delta2)^(E(i,3))<0 || 1-(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1
                            k=k+1;
                        end
                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*V(errenkada2)...
                            +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1)))*V(errenkada4));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*V(errenkada2)...
                            +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3)))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1)))*V(errenkada4);
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)<0 ||1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)>1
                            k=k+1;
                        end
                    end
                % i+k<N bada.
                else
                    % i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci1 *(E(i,3)/N)+lambda*(pi*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada2)...
                            +nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +(1-pi*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)*(1-delta2)^(E(i,3)))*V(errenkada4));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci1 *(E(i,3)/N)+lambda*(pi*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada2)...
                            +nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +(1-pi*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)*(1-delta2)^(E(i,3)))*V(errenkada4));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if (1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))<0 || (1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))>1 || 1-(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))<0 || 1-(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))>1 || nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))<0 || nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))>1
                            k=k+1;
                        end
                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada2)...
                            +nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)^(E(i,1))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)^(E(i,1))*(1-delta2)^(E(i,3)))*V(errenkada4));
                
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,3))*V(errenkada2)...
                            +nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)^(E(i,1))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +(1-nchoosek(E(i,1),1)*pi*(1-pi)^(E(i,1)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)^(E(i,1))*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-pi)^(E(i,1))*(1-delta2)^(E(i,3)))*V(errenkada4));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))>1 ||nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))<0 || nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-alpha)^(N-E(i,1)-E(i,3))-nchoosek(N-E(i,1)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,3)-1)*(1-delta2)^(E(i,3))>1
                            k=k+1;
                        end
                    end
                end
                if V0<V1
                    Vberria(i)=V0;
                    politika(i)=0;
                else
                    Vberria(i)=V1;
                    politika(i)=1;
                end
            end

            % (0,j,k) egoera non j+k<=N eta j,k>=1.
            if E(i,1)==0 && E(i,2)~=0 && E(i,3)~=0
                [~, errenkada1] = ismember([0,E(i,2)-1,E(i,3)], E, 'rows');
                [~, errenkada2] = ismember([0,E(i,2),E(i,3)-1], E, 'rows');
                [~, errenkada3] = ismember([1,E(i,2),E(i,3)], E, 'rows');
                [~, errenkada4] = ismember([0,E(i,2),E(i,3)], E, 'rows');
                politika(i)=0;
                % j+k=N bada.
                if E(i,2)+E(i,3)==N
                    Vberria(i)= ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*V(errenkada2)...
                            +(1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2)))*V(errenkada4));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))>1 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))<0 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))>1
                        k=k+1;
                    end
                %j+k<N bada.
                else
                    Vberria(i)= ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))*V(errenkada2)...
                            +nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +(1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3)))*V(errenkada4));
                    % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                    if nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))>1 || nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,2)-E(i,3))-nchoosek(N-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,2)-E(i,3)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1
                        k=k+1;
                    end
                end
                if delta1*(1-delta2)*(1-alpha) <0 || delta1*(1-delta2)*(1-alpha)>1 || delta2*(1-delta1)*(1-alpha) <0 || delta2*(1-delta1)*(1-alpha)>1 || alpha*(1-delta1)*(1-delta2)<0|| alpha*(1-delta1)*(1-delta2)>1|| 1-delta1*(1-delta2)*(1-alpha)-delta2*(1-delta1)*(1-alpha)-alpha*(1-delta1)*(1-delta2) <0 || 1-delta1*(1-delta2)*(1-alpha)-delta2*(1-delta1)*(1-alpha)-alpha*(1-delta1)*(1-delta2)>1 
                    k=k+1;
                end
            end

            % (i,j,k) egoera non i+j+k<=N eta i,j,k>=1.
            if E(i,1)~=0 && E(i,2)~=0 && E(i,3)~=0
                [~, errenkada1] = ismember([E(i,1)-1,E(i,2),E(i,3)], E, 'rows');
                [~, errenkada2] = ismember([E(i,1)-1,E(i,2)+1,E(i,3)], E, 'rows');
                [~, errenkada3] = ismember([E(i,1)-1,E(i,2),E(i,3)+1], E, 'rows');
                [~, errenkada4] = ismember([E(i,1),E(i,2)-1,E(i,3)], E, 'rows');
                [~, errenkada5] = ismember([E(i,1),E(i,2),E(i,3)-1], E, 'rows');
                [~, errenkada6] = ismember([E(i,1)+1,E(i,2),E(i,3)], E, 'rows');
                [~, errenkada7] = ismember([E(i,1),E(i,2),E(i,3)], E, 'rows');
                % i+j+k=N bada.
                if E(i,1)+E(i,2)+E(i,3)==N
                    % i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2)))*V(errenkada7));

                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2)))*V(errenkada7));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if (1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || (1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || 1-(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || 1-(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))>1 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))<0 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))>1
                            k=k+1;
                        end

                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2)))*V(errenkada7));

                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*V(errenkada5)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2)))*V(errenkada7));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))>1
                            k=k+1;
                        end
                    
                    end
                
                % i+j+k<N bada.
                else
                    % i=1 bada.
                    if E(i,1)==1
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada5)...
                            +nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada6)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3)))*V(errenkada7));
                   
                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada1)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada2)...
                            +(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada5)...
                            +nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada6)...
                            +(1-pi*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-pi)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3)))*V(errenkada7));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if (1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || (1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || 1-(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || 1-(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || 1-gamma*(E(i,2)/(N-1))*rho*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-gamma*(E(i,2)/(N-1))*(1-rho)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*(1-gamma*(E(i,2)/(N-1)))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1
                            k=k+1;
                        end
                    % i>1 bada.
                    else
                        pi=0;
                        V0= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada5)...
                            +nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada6)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3)))*V(errenkada7));

                        pi=1;
                        V1= cv*(E(i,1)/N)*pi+ci0*(E(i,2)/N)+ci1*(E(i,3)/N)+lambda*(nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada1)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada2)...
                            +nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada3)...
                            +nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada4)...
                            +nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))*V(errenkada5)...
                            +nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*V(errenkada6)...
                            +(1-nchoosek(E(i,1),1)*pi*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*rho*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*(1-pi)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-pi)*(1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3)))*V(errenkada7));
                        % Probabilitateak ondo definituta daudela egiaztatuko dugu.
                        if nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0 || nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))<0  || nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))>1 || nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0  || nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))<0 || 1-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*rho*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,1),1)*gamma*(E(i,2)/(N-1))*(1-rho)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1)-1)*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,2),1)*delta1*(1-delta1)^(E(i,2)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta2)^(E(i,3))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(E(i,3),1)*delta2*(1-delta2)^(E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3))-nchoosek(N-E(i,1)-E(i,2)-E(i,3),1)*alpha*(1-alpha)^(N-E(i,1)-E(i,2)-E(i,3)-1)*((1-gamma*(E(i,2)/(N-1))))^(E(i,1))*(1-delta1)^(E(i,2))*(1-delta2)^(E(i,3))>1
                            k=k+1;
                        end
                    end
                end
                if V0<V1
                    Vberria(i)=V0;
                    politika(i)=0;
                else
                    Vberria(i)=V1;
                    politika(i)=1;
                end
            end
        end
end
% Probabilitateak ez badaude ondo definituta mezu bat erakutsi.
if k>0
    disp("Probabilitateak ez daude ondo definituta daude!")
end

%Egoera bakoitzaren alboan egoera horren balio optimoa.
[E,Vberria]

% Egoera bakoitzean politika optimoaren grafikoa egin.
figure;
koloreak = [0 0 1; 0 0 0.1]; % Puntuen koloreak: urdina politika 0 eta beltza politika 1.
scatter3(E(:,1)/N,E(:,2)/N,E(:,3)/N, 20, koloreak(politika + 1, :),'filled');
%title("Politika optimoa $S$, $I_0$ eta $I_1$-ren proportzioen arabera","Interpreter","latex");
xlabel("$m_S(t)$","Interpreter","latex");
ylabel("$m_{I_0}(t)$","Interpreter","latex");
zlabel("$m_{I_1}(t)$","Interpreter","latex");
hold on; 
politika0 = scatter(nan, nan, 20, koloreak(1, :), 'filled');
politika1= scatter(nan, nan, 20, koloreak(2, :), 'filled');
hold off;
legend([politika0, politika1], {'Politika optimoa 0', 'Politika optimoa 1'}, "Interpreter", "latex"); 
grid on;

end