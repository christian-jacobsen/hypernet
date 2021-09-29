import numpy as np
from scipy.integrate import ode

from const import const


class ODESolver(object):

    ###########################################################################
    ### Initialization
    ###########################################################################
    def __init__(self, processes):

        # Initial/Final time
        self.t0     = 0.
        self.t_end  = 1.

        # Integrator time step
        self.dt0    = 1.e-14
        self.dt_str = 1.015
        self.dt_max = 1.e-4

        # Reactive processes
        self.processes = processes


    ###########################################################################
    ### Solve ODE
    ###########################################################################
    def f(self, t, y, arg):
        '''
        Input matrices:
        - arg[0] = Excitation rates
        - arg[1] = Dissociation rates
        - arg[2] = Recombination rates
        '''
        exc_term = np.matmul(arg[0], y) * y[-1]
        dis_term = np.matmul(arg[1], y) * y[-1]
        rec_term = np.squeeze(arg[2]) * np.power(y[-1], 3)
        return exc_term + dis_term + rec_term

    def jac(self, t, y, arg):
        '''Jacobian calculation: jac[i,j] = df[i] / dy[j].'''
        d = np.shape(y)[0]
        J = np.zeros((d,d), dtype=np.float64)
        J[:,:-1] = (arg[0][:,:-1] + arg[1][:,:-1]) * y[-1]
        J[:, -1] = np.matmul(arg[0][:,:-1], y[:-1]) + np.matmul(arg[1][:,:-1], y[:-1]) \
                    + 3*np.squeeze(arg[2])*np.power(y[-1], 2)
        return J

    def eval_J(self, t, y, arg):
        '''Eigenvalues of the Jacobian.'''
        J = self.jac(t, y, arg)
        # norm = np.hstack(([np.linalg.norm(i) for i in np.split(J, J.shape[1], axis=1)]))
        norm = np.linalg.norm(J)
        eigen, _ = np.linalg.eig(J)
        return norm, eigen

    def solve(self, rho_0, params, eval_jac=False, eval_f=False):
        '''Solving the ODE.'''

        r = ode(self.f, self.jac).set_integrator('vode', method='bdf', with_jacobian=True, atol=1.e-20)
        r.set_initial_value(rho_0, self.t0).set_f_params(params).set_jac_params(params)

        # Appending Data
        t = np.array([[self.t0]])
        y = np.expand_dims(rho_0,0)
        if eval_jac:
            norm_J, eig_J = self.eval_J(self.t0, rho_0, params)
        if eval_f:
            f_rhs = self.f(self.t0, rho_0, params)

        dt = self.dt0
        while r.successful() and r.t <= self.t_end:
            r.integrate(r.t+dt)
            if eval_jac:
                norm_J, eig_J = list(map(lambda x,y: np.vstack((x,y)), \
                    [norm_J, eig_J], self.eval_J(r.t, r.y, params)))
            if eval_f:
                f_rhs = np.vstack((f_rhs, self.f(r.t, r.y, params)))
            t  = np.vstack((t, np.expand_dims(r.t,0)))
            y  = np.vstack((y, np.expand_dims(r.y,0)))
            dt = min(dt*self.dt_str, self.dt_max)
        assert (y > 0.).all()

        out = [t, y]
        if eval_jac:
            out.extend([norm_J, eig_J])
        if eval_f:
            out.append(f_rhs)

        return out

    def equilibrium(self, rho):
        '''ODE Final Equilibrium Condition.'''

        # Translatiional partition functions
        Qtr_O  = np.power( 2*np.pi*const.UKb*T*const.mO / const.hp**2, 3./2.)
        Qtr_O2 = np.power( 2., 3./2.) * Qtr_O

        # Internal partition functions
        Q_O    = 9
        Q_O2   = Q_bins

        # Calculate rho_O
        mass = np.sum(rho)
        C1 = 2 / const.mO * Qtr_O2 / (Q_O * Qtr_O)**2
        C2 = 1 / np.sum(Q_O2)
        C3 = - mass / np.sum(Q_O2)
        
        r = np.roots([C1, C2, C3])
        rho_O = r[r>0]

        # Calculate rho_O2 bins
        rho_O2 = (mass - rho_O) * Q_O2 / np.sum(Q_O2)

        return np.append(rho_O2, rho_O)


    ###########################################################################
    ### Master Equation matrices
    ###########################################################################
    def get_matrices(self, rates_coeff, Q_bins, T):
        '''Get all the rates matrices for the ODE.'''
        K_e = self.excit(self.processes['excit'], rates_coeff, Q_bins) / const.mO
        K_d = self.disso(self.processes['diss'], rates_coeff) / const.mO
        K_r = self.recom(self.processes['recomb'], rates_coeff, Q_bins, T) / const.mO**2 * 2
        return K_e, K_d, K_r

    def excit(self, mask, rates_coeff, Q_bins):
        '''Excit. & Relax. rates matrix'''

        # Obtain relaxation rates
        rates = rates_coeff['exch'] + rates_coeff['inel']
        idx   = list(map(tuple, np.transpose(np.where(rates == 0.0))))
        for i,j in idx:
            if i != j:
                rates[i,j] = rates[j,i] * Q_bins[j] / Q_bins[i]

        # Construct Excit. & Relax. matrix
        d = np.shape(rates)[0] + 1
        K = np.zeros((d,d), dtype=np.float64)
        if mask:
            K[:-1,:-1] = -np.diag(np.sum(rates, axis=1)) + np.transpose(rates)
        return K

    def disso(self, mask, rates_coeff):
        '''Dissociation rates matrix'''

        # Obtain Dissociation rates
        rates = rates_coeff['diss']

        # Construct Dissociation matrix
        d = np.shape(rates)[0] + 1
        K = np.zeros((d,d), dtype=np.float64)
        if mask:
            K[:-1,:-1] = np.diag(-rates)
            K[ -1,:-1] = rates
        return K

    def recom(self, mask, rates_coeff, Q_bins, T):
        '''Recombination rates matrix'''

        # Obtain Recombination rates
        rates = np.multiply(rates_coeff['diss'], self.Keq_diss(Q_bins, T))

        # Construct Recombination matrix
        d = np.shape(rates)[0] + 1
        K = np.zeros((d,1), dtype=np.float64)
        if mask:
            K[:-1,0] = rates
            K[ -1,0] = -np.sum(rates)
        return K

    def Keq_diss(self, Q_bins, T):
        '''Diss.-Recomb. equilibrium constants matrix'''

        # Translational partition functions
        Qtr_O  = np.power( 2*np.pi*const.UKb*T*const.mO / const.hp**2, 3./2.)
        Qtr_O2 = np.power( 2., 3./2.) * Qtr_O

        # Internal partition functions
        Q_O    = 9.
        Q_O2   = Q_bins

        return Q_O2 * Qtr_O2 / (Q_O * Qtr_O)**2
