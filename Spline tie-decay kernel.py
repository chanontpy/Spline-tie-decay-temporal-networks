import sympy as sp

ta, h, e1, e2, a, k = sp.symbols('t_a h epsilon_1 epsilon_2 a k')
vector = sp.Matrix([a,e1,k,e2])#at time ta, the weight takes value a
matrix = sp.Matrix([[ta**3, ta**2, ta, 1],
                    [3*(ta**2), 2*ta, 1, 0],
                    [(ta+h)**3,(ta+h)**2,ta+h,1],
                    [3*(ta+h)**2,2*(ta+h),1,0]])
coeff1 = sp.lambdify([a,h,k,e1,e2],matrix.solve(vector)[0])
coeff2 = sp.lambdify([a,h,ta,k,e1,e2],matrix.solve(vector)[1])
coeff3 = sp.lambdify([a,h,ta,k,e1,e2],matrix.solve(vector)[2])
coeff4 = sp.lambdify([a,h,ta,k,e1,e2],matrix.solve(vector)[3])

class CubicSpline:
    
    def __init__(self,init_val,step_size,init_time,fin_val,init_deriv,fin_deriv):
        self.init_val = init_val
        self.step_size = step_size
        self.init_time = init_time
        self.fin_val = fin_val
        self.init_deriv = init_deriv
        self.fin_deriv = fin_deriv
    
    def polynomial(self,t):
        a = self.init_val
        h = self.step_size
        ta = self.init_time
        k = self.fin_val
        e1 = self.init_deriv
        e2 = self.fin_deriv
        return coeff1(a,h,k,e1,e2)*(t**3) + coeff2(a,h,ta,k,e1,e2)*(t**2) + coeff3(a,h,ta,k,e1,e2)*t + coeff4(a,h,ta,k,e1,e2)

class ExponentialDecay(CubicSpline):
    
    def __init__(self,init_val,step_size,init_time,fin_val,init_deriv,fin_deriv,alpha,next_in_time):
        super().__init__(init_val,step_size,init_time,fin_val,init_deriv,fin_deriv)
        self.alpha = alpha
        self.next_in_time = next_in_time
        a = self.init_val
        h = self.step_size
        ta = self.init_time
        k = self.fin_val
        e1 = self.init_deriv
        e2 = self.fin_deriv
        self.recent_val = CubicSpline(a,h,ta,k,e1,e2).polynomial(next_in_time)
        
    def decay(self,t):
        return self.recent_val*sp.exp(self.alpha*self.next_in_time)*sp.exp(-self.alpha*t)
