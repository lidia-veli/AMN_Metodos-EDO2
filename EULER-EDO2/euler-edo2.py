import numpy as np

def euler_edo2(Fx, a, b, y0, y1, N, Px, Qx):
    '''funcion que resuleve la EDO2 de la forma y'' + p(x)y' + q(x)y = f(x)
    
    INPUTS:
    - Fx: funcion f(x)
    - a, b: extremos del intervalo 
    - y0, y1: condiciones iniciales para y(a) y y'(a)
    - N: numero de subintervalos
    - Px, Qx: funciones p(x) y q(x)
    '''

    h = (b-a)/N  # paso

    x = [a]
    u = [y0]  # u(x) = y(x)
    v = [y1]  # v(x) = y'(x)

    for k in range(N):
        x.append(x[k] + h)
        u.append(u[k] + h*v[k])
        v.append(v[k] + h*( - Px(x[k])*v[k] - Qx(x[k])*u[k] + Fx(x[k]))  )
    
    return x, u, v
