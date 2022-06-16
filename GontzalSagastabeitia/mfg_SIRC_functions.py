class Parameters:
    
    def __init__(self,N=float('inf')):
        """
        Default values
        """
        if N<float('inf'):
            self.gamma = 0.9     # infection rate
            self.rho = 0.8    # recovery rate
            self.SIR0 = (round(N/2),N-round(N/2),0)  #initial state of populations
            self.T = 3           # time horizon
            self.epsilon = 0          # value of ε
            self.ci = 500      # cost per unit time of infection
            self.N = N
            self.C = 10*N
        else:
            self.gamma = 0.9     # infection rate
            self.rho = 0.3    # recovery rate
            self.SIR0 = (3/9,1-3/9,0)  #initial state of populations (mI decreasing)
            self.T = 100           # time horizon
            self.epsilon = 0          # value of ε
            self.ci = 1000      # cost per unit time of infection
            self.A = [0,1]     #possible actions of player 0

#N=Inf case functions ------------------------------------------------------------------------------------------------

def mSIR(parameters,pi):
    if parameters.T==0:
        return parameters.SIR0
    else:
        mS=[parameters.SIR0[0]]
        mI=[parameters.SIR0[1]]
        mR=[parameters.SIR0[2]]
        for i in range(1,parameters.T+1):
            Si=mS[i-1]-parameters.gamma*mS[i-1]*mI[i-1]*pi[i-1]
            Ii=mI[i-1]+parameters.gamma*mS[i-1]*mI[i-1]*pi[i-1]-parameters.rho*mI[i-1]
            mS.append(Si)
            mI.append(Ii)
            mR.append(1-Si-Ii) #The proportion of R is always the remaining population
        return (mS,mI,mR)

def getVI(parameters):
    VI=[0]
    for i in range(0,parameters.T):
        VI.append(VI[i]+parameters.ci*(1-parameters.rho)**i)
    VI.reverse()
    return VI

def getVS(parameters,pi,BR=False):
    VS=[0]
    BR0=[]
    VI=getVI(parameters)
    mI=mSIR(parameters,pi)[1]
    for i in range(0,parameters.T):
        VaS=[]
        for a in parameters.A:
            VaS.append(f(a)+(1-parameters.gamma*mI[parameters.T-i-1]*a)*VS[0]+parameters.gamma*mI[parameters.T-i-1]*a*VI[parameters.T-i])
        VS=[min(VaS)]+VS
        if BR:
            BR0=[parameters.A[VaS.index(VS[0])]]+BR0
    if BR:
        return (BR0,VS)
    else:
        return VS

def f(a):
    return 2-a

def pi_cost(parameters,pi):
    """
    Computes the cost of π by adding the cost of every t=0,...,T
    """
    (mS,mI,mR)=mSIR(parameters,pi)
    cost=0
    for t in range(parameters.T+1):
        cost+=mS[t]*f(pi[t])+mI[t]*parameters.ci
    return cost

def compute_MFE_cost(parameters):
    """
    Computes the cost of the MFE. For that it computes MFE starting from π,
    then it uses the function pi_cost() to compute the cost.
    """
    i=0
    pi=pi_i=[1]*(parameters.T)
    pi0=getVS(parameters,pi,True)[0]
    while pi0!=pi:
        pi=pi0
        i+=1
        pi0=getVS(parameters,pi,True)[0]
        if i==100:
            break
    iterazioak=i if i<100 else print('Warning: MFE not found')
    mfe_cost=pi_cost(parameters,pi0+[0])
    return (mfe_cost,pi0)

def compute_OPT_cost(parameters):
    """
    Computes the cost of the OPT population strategy and the strategy itself.
    It uses the function pi_cost().
    We are assuming the optimal π has at most one jump.
    """
    pi_i=[1]*(parameters.T+1)
    pi_OPT=[pi_i]
    OPT_cost=pi_cost(parameters,pi_i)
    for i in range(1,parameters.T+2):
        pi_i=[0]*i+[1]*(parameters.T+1-i)
        cost_i=pi_cost(parameters,pi_i)
        if cost_i<OPT_cost:
            pi_OPT=[pi_i]
            OPT_cost=cost_i
        elif cost_i==OPT_cost:
            pi_OPT.append(pi_i)
    #print(str(len(pi_OPT))+' optimal strategies were found.')
    return (OPT_cost,pi_OPT)



#N<Inf case functions ------------------------------------------------------------------------------------------------
import numpy as np

def comp_prob(parameters,Ms,Mi,a,b):
    """
    Computes the transition probabilities at a given time t with populations Ms & Mi and policies BR(t)=a & π(t)=b
    """
    #player 0 probabilities:
    px=np.zeros(2)
    pm=np.zeros(3)
    
    C=parameters.C
    gamma=parameters.gamma/C
    rho=parameters.rho/C
    
    px=[gamma*Mi*a/(parameters.N), #prob of infection
        rho] #prob of recovery
    
    #mass probabilities:
    pm=[gamma*Ms*Mi*b/(parameters.N**2), #prob of inf when player 0 is NOT inf
        rho*Mi/(parameters.N), #prob of recovery
        gamma*Ms*(Mi+1)*b/(parameters.N**2)]#prob of inf when player 0 is inf
    return (px,pm)

def BR_N_players(parameters,pi):
    """
    Computes BR(pi), Vs^{BR,pi} and Vi^{BR,pi}
    """
    C=parameters.C
    Vs=np.zeros((C*parameters.T+1,parameters.N+1,parameters.N+1))
    Vi=np.zeros((C*parameters.T+1,parameters.N+1,parameters.N+1))
    BR=np.zeros((C*parameters.T+1,parameters.N+1,parameters.N+1))
    for t in range(C*parameters.T-1,-1,-1):
        for Ms in range(parameters.N+1):
            for Mi in range(parameters.N-Ms+1):
                if parameters.gamma*Mi*(Vi[t+1,Ms,Mi]-Vs[t+1,Ms,Mi])/parameters.N < 1:
                    BR[t,Ms,Mi]=1/C
                else:
                    BR[t,Ms,Mi]=parameters.epsilon
                px,pm=comp_prob(parameters,Ms,Mi,BR[t,Ms,Mi],pi[t,Ms,Mi])
                if (Ms,Mi)==(0,0):
                    Vs[t,Ms,Mi]=0
                    Vi[t,Ms,Mi]=0
                else:
                    Vs[t,Ms,Mi]=2-BR[t,Ms,Mi]+\
                                  px[0]*Vi[t+1,Ms,Mi]+\
                                  (pm[1]*Vs[t+1,Ms,Mi-1] if Mi>0 else 0) +\
                                  (pm[0]*Vs[t+1,Ms-1,Mi+1] if Ms>0 else 0) +\
                                 ((1-px[0]-pm[1]-pm[0])*Vs[t+1,Ms,Mi])
                    Vi[t,Ms,Mi]=parameters.ci+\
                        (pm[1]*Vi[t+1,Ms,Mi-1] if Mi>0 else 0)+\
                        (pm[2]*Vi[t+1,Ms-1,Mi+1] if Ms>0 else 0)+\
                        (1-px[1]-pm[1]-pm[2])*Vi[t+1,Ms,Mi]

    return (BR,Vs,Vi)

def compute_MFE_N_players(parameters,max_iterations):
    """
    Tries to find a fixed point of BR. If no fixed point found after max_iterations it stops.
    """
    N = parameters.N
    T = parameters.T
    C = parameters.C
    # The initial policy of the mass is to always confine.
    pi = np.zeros((C*parameters.T+1,parameters.N+1,parameters.N+1))

    BR, Vs, Vi = BR_N_players(parameters,pi)

    iteration_numbers = 1
    while iteration_numbers<max_iterations:
        number_of_differences = np.sum(1-np.isclose(BR, pi))
        if number_of_differences == 0:
            print("Equilibrium found after "+str(iteration_numbers)+" iterations.")
            return BR, Vs, Vi
        else:
            if np.random.rand() >= 0.1:
                pi = BR    
            else:
                pi = BR*0.1+pi*0.8
            BR, Vs, Vi = BR_N_players(parameters,pi)
            iteration_numbers += 1
    print("Equilibrium NOT found after "+str(iteration_numbers)+" iterations.")

def BR_cost_N_players(parameters,pi,Vs=None,Vi=None):
    """
    Returns the cost of BR(π). If no Vs or Vi are given as parameters, they are computed inside de function.
    """
    if Vs is None or Vi is None:
        BR,Vs,Vi=BR_N_players(parameters,pi)
    (Ms0,Mi0,Mr0)=parameters.SIR0
    N=parameters.N
    cost=Vs[0,Ms0,Mi0]*(Ms0/N)+Vi[0,Ms0,Mi0]*(Mi0/N)
    return cost



