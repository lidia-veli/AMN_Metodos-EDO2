def runge_kutta_edo2(Fx, a, b, y0, y1, N, Px, Qx):
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
        k11 = h*v[k]
        k12 = h*Fx(x[k], u[k], v[k])
        k21 = h*(v[k] + k12/2)
        k22 = h*Fx(x[k] + h/2, u[k] + k11/2, v[k] + k12/2)
        k31 = h*(v[k] + k22/2)
        k32 = h*Fx(x[k] + h/2, u[k] + k21/2, v[k] + k22/2)
        k41 = h*(v[k] + k32)
        k42 = h*Fx(x[k] + h, u[k] + k31, v[k] + k32)

        # actualizamos los vectores
        u.append(u[k] + (k11 + 2*k21 + 2*k31 + k41)/6)
        v.append(v[k] + (k12 + 2*k22 + 2*k32 + k42)/6)
        x.append(x[k] + h)
    
    return x, u, v



def rg_edo2(f, x, u, v, h, n):
    r = []
    t = []
    for i in range(n):

        k11=v
        k12= f(x,u,v)
        k21= v + (h/2)*k12
        k22= f(x+(h/2), u+(h/2)*k11, v+(h/2)*k12)
        k31= v + (h/2)*k22
        k32= f(x+(h/2), u+(h/2)*k21, v+(h/2)*k22)
        k41= v + h*k32
        k42= f(x+h, u+h*k31, v+h*k32)

        u = u + h * (((1/6) * k11) + ((1/3) * k21) + ((1/3) * k31) + ((1/6) * k41))
        v= v + h * (((1/6) * k12) + ((1/3) * k22) + ((1/3) * k32) + ((1/6) * k42))
        x = x + h
        r.append(x)
        t.append(u)
    return r, t
