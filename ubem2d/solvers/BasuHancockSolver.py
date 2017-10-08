import numpy as np
import numpy.linalg as nla
import scipy.linalg as sla
from ubem2d.fluids.BasicFlows import velocity_uniform_flow
from ubem2d.panel.PanelInfluence import velocity_source_panel
from ubem2d.panel.PanelInfluence import velocity_vortex_panel
from ubem2d.panel.BodyInfluence import source_influence_matrices_body
from ubem2d.panel.BodyInfluence import vortex_influence_matrices_body
from ubem2d.panel.BodyInfluence import velocity_source_body
from ubem2d.panel.BodyInfluence import velocity_vortex_body
from ubem2d.Errors import SolverError
from .HessSmithSolver import solve_hess_smith_body

__all__ = ['BasuHancockSolver']

class BasuHancockSolver():
    '''
    This class implements the unsteady boundary-element method, described in
    Basu and Hancock (JFM 1978), for flow past an airfoil.
    '''
    def __init__(self, body, wake, xref = -10, yref = 0, nref = 20,
        maxiters = 200, tol = 1.e-6, maxerr = 1.e-5, wakep_free = True,
        wake_body = True, wake_self = True):

        super().__init__()

        self._steps = 0                 # Number of steps taken thus far

        self._body = body               # Airfoil object
        self._wake = wake               # Wake object
        self._xref = xref               # x coordinate of potential ref point
        self._yref = yref               # y coordinate of potential ref point
        self._nref = nref               # Number of panels from ref point to LE
        self._maxiters = maxiters       # Max number of wake panel iterations
        self._tol = tol                 # Unsteady convergence tolerance
        self._maxerr = maxerr           # Max allowed Neumann/Kutta error
        self._wakep_free = wakep_free   # False: wake panel bisects TE
        self._wake_body = wake_body     # Whether body influences wake
        self._wake_self = wake_self     # Whether wake influences itself

        self._delk = None               # Wake panel length
        self._thk = None                # Wake panel inclination to +x-axis
        self._shed_circ = None          # Circulation of last shed vortex
        self._shed_x = None             # x coordinate of last shed vortex
        self._shed_y = None             # y coordinate of last shed vortex

        # Construct body influence matrices and perform LU factorization
        self._At, self._An = source_influence_matrices_body(self._body)
        self._Bt, self._Bn = vortex_influence_matrices_body(self._body)
        self._lup = sla.lu_factor(self._An)

    def step(self, dt = 0, uinf=(1,0)):
        if (self._steps == 0):
            return self.steady_step(uinf)
        elif (dt != 0):
            return self.unsteady_step(dt, uinf)
        else:
            raise ValueError('Unsteady time step must be nonzero')

    def steady_step(self, uinf):
        soln = solve_hess_smith_body(uinf, self._body, self._At, self._An,
            self._Bt, self._Bn)
        phik = self.compute_potential(uinf, soln.qt, soln.sigma, soln.gamma, 0)
        # Make initial guess for wake panel length and inclination
        self._delk = self._body.perimeter/self._body.nedge
        self._thk = self.trailing_edge_bisector()
        return self.post_step(soln.sigma, soln.gamma, phik, soln.cp, 0, 0, 0)

    def unsteady_step(self, dt, uinf):
        # Kinematic update and Kutta condition
        dxdt = (self._body._xmid - self._xmid)/dt
        dydt = (self._body._ymid - self._ymid)/dt
        vn = dxdt*self._body.nx + dydt*self._body.ny
        sigk,gamk,uwk,vwk = self.solve_implicit_kutta(uinf,vn,dt)
        # Kelvin circulation theorem
        L = self._body.perimeter
        gamwk = (L/self._delk)*(self._gam-gamk)
        # Check Neumann boundary condition and Kutta condition
        qt,qn = self.flow(uinf, sigk, gamk, gamwk)
        q = np.sqrt(qt*qt + qn*qn)
        error_neumann = nla.norm(qn - vn)
        error_kutta = np.abs(q[0]**2 - q[-1]**2 - 2*L*(gamk-self._gam)/dt)
        if (error_neumann > self._maxerr):
            raise SolverError('Neumann error: {}'.format(
                error_neumann))
        if (error_kutta > self._maxerr):
            raise SolverError('Kutta error: {}'.format(
                error_kutta))
        # Compute potential and pressure distribution via unsteady Bernoulli
        phik = self.compute_potential(uinf, qt, sigk, gamk, gamwk)
        dphidt = (phik-self._phi)/dt
        spdinf = nla.norm(uinf)
        cp = 1.-(q*q + 2*dphidt)/spdinf**2
        # Detach wake panel and advect the wake
        shed_circ = gamwk*self._delk
        shed_x = self._body.x[0] + .5*self._delk*np.cos(self._thk) + uwk*dt
        shed_y = self._body.y[0] + .5*self._delk*np.sin(self._thk) + vwk*dt
        self._wake.append(shed_circ, shed_x, shed_y)
        self.advect_wake(uinf, sigk, gamk, dt)
        return self.post_step(sigk, gamk, phik, cp, shed_circ, shed_x, shed_y)

    def post_step(self, sigk, gamk, phik, cp, shed_circ, shed_x, shed_y):
        '''
        Update values and store for use in the next time step.  Return the 
        solution data.
        '''
        self._steps += 1
        self._circ_bound = gamk*self._body._perimeter
        self._xmid = self._body.xmid
        self._ymid = self._body.ymid
        self._sig = sigk
        self._gam = gamk
        self._phi = phik
        self._cp = cp
        return (self._sig, self._gam, self._cp, shed_circ, shed_x, shed_y)

    def trailing_edge_bisector(self):
        '''
        Return the angle in (-pi,pi) of an aft-pointing unit vector which
        bisects the trailing-edge angle.
        '''
        dx = .5*(self._body.tx[-1] - self._body.tx[0])
        dy = .5*(self._body.ty[-1] - self._body.ty[0])
        return np.arctan2(dy,dx)

    def solve_implicit_kutta(self, uinf, vn, dt):
        (uinft,uinfn) = self.flow_onset(uinf)
        (Wvt,Wvn) = self.flow_wake()
        L = self._body.perimeter
        n = self._body.nedge
        if (not self._wakep_free):
            self._thk = self.trailing_edge_bisector()

        # Wake panel iteration
        converged = False
        for i in range(self._maxiters):
            (Wpt,Wpn) = self.flow_wake_panel(1.)
            bk = (L/self._delk)*Wpn - np.sum(self._Bn,1)
            ck = -uinfn - (L/self._delk)*self._gam*Wpn - Wvn + vn
            xxk = sla.lu_solve(self._lup, bk)
            yyk = sla.lu_solve(self._lup, ck)
            alpha1 = np.dot(self._At[0,:],xxk) + np.sum(self._Bt[0,:]) \
                - (L/self._delk)*Wpt[0]
            beta1 = np.dot(self._At[0,:],yyk) \
                + (L/self._delk)*self._gam*Wpt[0] + Wvt[0] + uinft[0]
            alphaN = np.dot(self._At[-1,:],xxk) + np.sum(self._Bt[-1,:]) \
                - (L/self._delk)*Wpt[-1]
            betaN = np.dot(self._At[-1,:],yyk) \
                + (L/self._delk)*self._gam*Wpt[-1] + Wvt[-1] + uinft[-1]
            zeta = alpha1**2 - alphaN**2
            eta = 2*(alpha1*beta1 - alphaN*betaN - L/dt)
            chi = beta1**2 - betaN**2 + 2*L*self._gam/dt + (vn[0])**2 \
                - (vn[-1])**2

            # Solve quadratic equation for gamk
            gamk_vals = np.roots((zeta, eta, chi))
            # Choose root with the smallest absolute value
            if np.abs(gamk_vals[0]) < np.abs(gamk_vals[1]):
                gamk = gamk_vals[0]
            else:
                gamk = gamk_vals[1]
            sigk = gamk*xxk + yyk

            # Compute resulting flow at wake panel midpoint
            xwkmid = np.array([self._body.x[0] + \
                .5*self._delk*np.cos(self._thk)])
            ywkmid = np.array([self._body.y[0] + \
                .5*self._delk*np.sin(self._thk)])
            (uv,vv) = velocity_vortex_panel(self._body.x[:-1],
                self._body.y[:-1], self._body.tx, self._body.ty,
                self._body.edge, gamk*np.ones(n), xwkmid, ywkmid)
            (us,vs) = velocity_source_panel(self._body.x[:-1],
                self._body.y[:-1], self._body.tx, self._body.ty,
                self._body.edge, sigk, xwkmid, ywkmid)
            (uw,vw) = self._wake.velocity(xwkmid, ywkmid)
            uwk = uv + us + uw + uinf[0]
            vwk = vv + vs + vw + uinf[1]

            # Update wake panel geometry and check for convergence
            self._delk = np.sqrt(uwk**2 + vwk**2)[0]*dt
            if (self._wakep_free):
                self._thk = np.arctan2(vwk,uwk)[0]
            if (i > 0 and nla.norm([uwk-uwk0,vwk-vwk0]) < self._tol):
                converged = True
                break
            uwk0 = uwk
            vwk0 = vwk
        # End of wake panel iteration loop
        if (not converged):
            raise SolverError('Unsteady wake panel failed to converge')
        return sigk,gamk,uwk,vwk

    def flow(self, uinf, sigk, gamk, gamwk):
        '''
        Return the net tangential and normal components of the flow at the
        panel midpoints.
        '''
        (uinft,uinfn) = self.flow_onset(uinf)
        (Pant,Pann) = self.flow_body_panels(sigk,gamk)
        (Wpt,Wpn) = self.flow_wake_panel(gamwk)
        (Wt,Wn) = self.flow_wake()
        return (
            uinft + Pant + Wpt + Wt,
            uinfn + Pann + Wpn + Wn)

    def flow_onset(self, uinf):
        '''
        Return the tangential and normal components of the flow at the panel
        midpoints due to the onset flow.
        '''
        return (uinf[0]*self._body.tx + uinf[1]*self._body.ty,
            uinf[0]*self._body.nx + uinf[1]*self._body.ny)

    def flow_body_panels(self,sigk,gamk):
        '''
        Return the tangential and normal components of the flow at the panel
        midpoints due to the source and vortex distributions.
        '''
        bodyt = np.dot(self._At,sigk) + gamk*np.sum(self._Bt,1)
        bodyn = np.dot(self._An,sigk) + gamk*np.sum(self._Bn,1)
        return (bodyt,bodyn)

    def flow_wake_panel(self,gamwk):
        '''
        Return the tangential and normal components of the flow at the panel
        midpoints due to the wake panel.
        '''
        x1 = self._body.x[0]
        y1 = self._body.y[0]
        tx = np.cos(self._thk)
        ty = np.sin(self._thk)
        (u,v) = velocity_vortex_panel(x1,y1,tx,ty,self._delk,
            gamwk, self._body.xmid, self._body.ymid)
        return (u*self._body.tx + v*self._body.ty,
            u*self._body.nx + v*self._body.ny)

    def flow_wake(self):
        '''
        Return the tangential and normal components of the flow at the panel
        midpoints due to the wake vortices.
        '''
        (u,v) = self._wake.velocity(self._body.xmid, self._body.ymid)
        return (u*self._body.tx + v*self._body.ty,
            u*self._body.nx + v*self._body.ny)

    def compute_potential(self, uinf, qt, sigk, gamk, gamwk):
        '''
        Compute the unsteady potential at the panel midpoints.  Do so by
        performing a line integral of the velocity field from the reference
        point up to the leading edge, and then continuing around the upper
        and lower surfaces.
        '''
        n = self._body.nedge
        # Sequence of panels from reference point to airfoil's leading edge
        (xle,yle) = self._body.leading_edge
        xpp = np.linspace(self._xref, xle, self._nref+1)
        ypp = np.linspace(self._yref, yle, self._nref+1)
        # Contributions from onset flow, body source & vortex panels, wake
        (u,v) = velocity_uniform_flow(uinf,xpp[:-1],ypp[:-1])
        (us,vs) = velocity_source_body(self._body,sigk,xpp[:-1],ypp[:-1])
        (uv,vv) = velocity_vortex_body(self._body,gamk*np.ones(n),
            xpp[:-1],ypp[:-1])
        (uw,vw) = self._wake.velocity(xpp[:-1],ypp[:-1])
        # Net flow, except wake panel
        u,v = u+us+uv+uw, v+vs+vv+vw
        # Contribution from wake panel
        if (self._delk is not None and self._thk is not None):
            x1 = self._body.x[0]
            y1 = self._body.y[0]
            tx = np.cos(self._thk)
            ty = np.sin(self._thk)
            (uwp,vwp) = velocity_vortex_panel(x1,y1,tx,ty,self._delk,gamwk,
                xpp[:-1],ypp[:-1])
            u,v = u+uwp, v+vwp
        # Line integral of flow from reference point to the leading edge
        le = self._body.le   # Index of the airfoil's leading-edge corner
        phi = np.zeros(n+1)  # We'll compute the potential at each corner
        phi[le] = np.sum(u*np.diff(xpp) + v*np.diff(ypp))
        # Line integral, starting from LE and marching along upper surface
        for i in range(le-1,-1,-1):
            phi[i] = phi[i+1] - qt[i]*self._body.edge[i]
        # Line integral, starting from LE and marching along lower surface
        for i in range(le+1,n+1):
            phi[i] = phi[i-1] + qt[i-1]*self._body.edge[i-1]
        # Potential at midpoints is the average of the potential at corners
        return .5*(phi[:-1]+phi[1:])

    def advect_wake(self, uinf, sigk, gamk, dt):
        '''
        Advect the wake vortices with a simple one-step, explicit Euler
        integration scheme.
        '''
        npan = self._body.nedge
        nvort = len(self._wake)
        vx = uinf[0]*np.ones(nvort)
        vy = uinf[1]*np.ones(nvort)
        if (self._wake_body):
            (us,vs) = velocity_source_body(self._body, sigk, self._wake.x,
                self._wake.y)
            (uv,vv) = velocity_vortex_body(self._body, gamk*np.ones(npan),
                self._wake.x, self._wake.y)
            vx += us + uv
            vy += vs + vv
        if (self._wake_self):
            (us,vs) = self._wake.self_velocity()
            vx += us
            vy += vs
        self._wake.advect(vx,vy,dt)
