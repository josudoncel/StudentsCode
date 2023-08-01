import numpy as np
#def phi(x):
    #return 1+x**2

def phi(x):
    return np.exp(x)

#def phi(x):
    #return 3*x**10


from scipy.optimize import minimize
from scipy.optimize import LinearConstraint



def funcion_obj(x_vals,*args):
    "c: lista con los costes de la funcion objetivo"
    "phi: funcion de latencia del problema"
    "x: estrategia que queremos que minimize la funcion objetivo"
    "s_others: lista con la estrategia del otro jugador"
    c,phi,s_others = args
    sumandos=[]
    for j in range(len(x_vals)):
        y=0
        for i in range(len(s_others)):
            y += s_others[i][j]
        sumandos.append(c[j]*phi(x_vals[j] + y)*x_vals[j])
    return sum(sumandos)
    
def restriccion(x_vals, f_i):
    "f_i:flujo total que tiene el jugador del que queremos calcular la mejor respuesta"
    return sum(x_vals) - f_i


#def best_response_player_i(c, phi, f_i, s_others, s_inicial_i):
    #"s_inicial_i: estrategia inicial de nuestro jugador tal que que verifica las restricciones"
    #rest = {'type': 'eq', 'fun': restriccion, 'args': (f_i,)}
    #sol = minimize(funcion_obj, s_inicial_i, args=(c, phi, s_others), constraints=[rest], method='trust-constr')
    #return sol.x

def best_response_player_i(c, phi, f_i, s_others, s_inicial_i):
    "s_inicial_i: estrategia inicial de nuestro jugador tal que que verifica las restricciones"
    rest = {'type': 'eq', 'fun': restriccion, 'args': (f_i,)}
    bounds = [(0, np.inf) for _ in range(len(s_inicial_i))]
    sol = minimize(funcion_obj, s_inicial_i, args=(c, phi, s_others), constraints=[rest], bounds=bounds, method='trust-constr')
    return sol.x

import random

def estrategia_inicial(f_i, l):
    lista = [abs(random.uniform(0, l)) for i in range(l)]
    suma = sum(lista)
    lista = [int(round(f_i * x/suma)) for x in lista]
    while sum(lista) != f_i:
        if sum(lista) < f_i:
            i = random.randint(0, l-1)
            lista[i] += 1
        else:
            i = random.randint(0, l-1)
            lista[i] -= 1
    return lista



import numpy as np

import numpy as np

def converge_to_nash_equilibrium(c, J, f, phi, s, max_iter=100, tolerancia=1e-3):
    """
    Ejecuta el algoritmo de mejor respuesta hasta converger al equilibrio de Nash.
    c: lista con los costes de cada ruta
    J: numero de jugadores
    f: lista con los flujos de cada jugador
    phi: funcion de latencia
    s: perfil estratégico inicial
    max_iter: número máximo de iteraciones permitido
    tolerancia: tolerancia para comparar perfiles estratégicos
    """
    iteraciones=0
    while True:
        iteraciones = iteraciones + 1
        ant=s.copy()
        for j in range(J):
            estrategia_inicial=s.pop(j)
            s.insert(j, best_response_player_i(c,phi,f[j], s, estrategia_inicial))# calcular su mejor respuesta al perfil estratégico actual
        if iteraciones == max_iter:
            print(s)#si se supera el maximo de iteraciones, devolvemos error
            raise ValueError("No se pudo converger al equilibrio de Nash después de {} iteraciones".format(max_iter))
            break
        if np.allclose(ant, s, rtol=tolerancia, atol=tolerancia):   # si no hubo cambios, el perfil actual es un equilibrio de Nash
            return [list(tuple(s[i])) for i in range(J)]


costes = [10,15,13]
numero_jugadores = 2
flujos = [3,4]
perfil_inicial = [[1,1,1],[1.5,0.5,2]]



print(converge_to_nash_equilibrium(costes,numero_jugadores,flujos,phi,perfil_inicial))
print(best_response_player_i(costes,phi,3,[[1.4309417809760168, 1.2526898117476346, 1.3163684072763484]],[0.5,1.5,1]))
print(best_response_player_i(costes,phi,4,[[1.0835240577188947, 0.9324646707457083, 0.9840112715353969]],[1,2,1]))
#print(best_response_player_i(costes,phi,4,[[20.567838541400405, 15.242942025753035, 17.63558751699991, 16.390673133275683, 10.309050024539351, 11.763915980191669, 8.089992777839946], [4.113943055513463, 3.0476927645141827, 3.530139942292061, 3.2759034756669876, 2.0625577999033684, 2.349670338502458, 1.6200926236074749],[39.08066526221182, 28.966332817180252, 33.4934610191893, 31.15463463887828, 19.58480532560425, 22.36137791758067, 15.358723019355411]], [1,1,1,0,0,0,0]))
#print(best_response_player_i(costes,phi,190,[[20.567838541400405, 15.242942025753035, 17.63558751699991, 16.390673133275683, 10.309050024539351, 11.763915980191669, 8.089992777839946], [4.113943055513463, 3.0476927645141827, 3.530139942292061, 3.2759034756669876, 2.0625577999033684, 2.349670338502458, 1.6200926236074749], [0.8233950694005003, 0.6094671419933956, 0.7071941799692962, 0.6543811226067241, 0.4134849579728359, 0.4676751761410225, 0.32440235191622535]],[37, 22, 44, 20, 20, 5, 42]))

