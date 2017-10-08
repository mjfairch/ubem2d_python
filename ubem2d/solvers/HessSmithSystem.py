import numpy as np
import numpy.linalg as nla
import scipy.linalg as sla
from ubem2d.panel.PanelInfluence import velocity_source_panel
from ubem2d.panel.PanelInfluence import velocity_vortex_panel
from ubem2d.panel.BodyInfluence import velocity_source_body
from ubem2d.panel.BodyInfluence import velocity_vortex_body

__all__ = ['HessSmithSystem']

class HessSmithSystem:
    def __init__(self, bodies):
        if (type(bodies) not in [list,tuple]):
            bodies = [bodies]
        Nb = len(bodies)
        Ns = [len(body) for body in bodies]

        # Order panel data from the first body through the last body. Then
        # compute start/end indices for each body.
        x1 = np.concatenate([body.x[:-1] for body in bodies])
        y1 = np.concatenate([body.y[:-1] for body in bodies])
        tx = np.concatenate([body.tx for body in bodies])
        ty = np.concatenate([body.ty for body in bodies])
        nx = np.concatenate([body.nx for body in bodies])
        ny = np.concatenate([body.ny for body in bodies])
        edge = np.concatenate([body.edge for body in bodies])
        xmid = np.concatenate([body.xmid for body in bodies])
        ymid = np.concatenate([body.ymid for body in bodies])
        a = np.concatenate([[0], np.cumsum(Ns)[:-1]]) # start indices
        b = np.cumsum(Ns) - 1                         # end indices

        # Compute tangential/normal source/vortex influence matrices
        N = sum(Ns)  # Total number of panels across all bodies
        At, An = np.zeros((N, N)), np.zeros((N, N))  # Source influences
        Bt, Bn = np.zeros((N, N)), np.zeros((N, N))  # Vortex influences
        for j in range(N):
            (u,v) = velocity_source_panel(x1[j], y1[j], tx[j], ty[j], edge[j],
                1., xmid, ymid)
            At[:,j], An[:,j] = u*tx + v*ty, u*nx + v*ny
            (u,v) = velocity_vortex_panel(x1[j], y1[j], tx[j], ty[j], edge[j],
                1., xmid, ymid)
            Bt[:,j], Bn[:,j] = u*tx + v*ty, u*nx + v*ny
        for j in range(N):
            # Update panel self-influences (based on hand computation)
            At[j,j], Bn[j,j] = 0, 0
            An[j,j], Bt[j,j] = .5, .5

        # Compute Hess-Smith matrix
        A = np.zeros((N+Nb, N+Nb))
        A[:-Nb,:-Nb] = An
        for k in range(Nb):
            A[:-Nb,N+k] = np.sum(Bn[:,a[k]:b[k]+1],1)
            A[N+k,:-Nb] = At[a[k],:] + At[b[k],:]
            #print('ERROR! Lower-right off-diagonal entries incorrectly zero.')
            A[N+k,N+k] = Bt[a[k],:].sum() + Bt[b[k],:].sum()

        # Store data for later usage
        self._bodies = bodies
        self._Nb = Nb
        self._N = N
        self._At = At
        self._An = An
        self._Bt = Bt
        self._Bn = Bn
        self._tx = tx
        self._ty = ty
        self._nx = nx
        self._ny = ny
        self._edge = edge
        self._a = a
        self._b = b
        self._LU = sla.lu_factor(A)

    def solve(self, uinf):
        '''
        Return source strengths along each body and circulation per unit
        length along each body.
        '''
        # Convenience variables
        N, Nb = self._N, self._Nb
        tx, ty, nx, ny = self._tx, self._ty, self._nx, self._ny
        a, b = self._a, self._b
        # Construct right-hand side
        rhs = np.zeros(N + Nb)
        rhs[0:N] = -(uinf[0]*nx + uinf[1]*ny)
        for k in range(Nb):
            rhs[N+k] = -(uinf[0]*(tx[a[k]] + tx[b[k]]) +
                uinf[1]*(ty[a[k]] + ty[b[k]]))
        # Solve the system
        soln = sla.lu_solve(self._LU, rhs)
        # Extract source and circulation strenghts for each body
        sigma = [soln[a[k]:b[k]+1] for k in range(Nb)]
        gamma = soln[N:]
        return (sigma, gamma)

    def flow_self(self, uinf, soln):
        '''
        Compute tangential flow at panel midpoints (normal flow is zero).
        '''
        At, Bt, tx, ty = self._At, self._Bt, self._tx, self._ty
        a, b = self._a, self._b
        (sigma, gamma) = soln
        qt = uinf[0]*tx + uinf[1]*ty  # due to onset flow
        qt += np.dot(At,np.concatenate(sigma))  # due to source terms
        for k in range(self._Nb):  # due to circulation round each body
            qt += gamma[k]*np.sum(Bt[:,a[k]:b[k]+1],1)
        return [qt[a[k]:b[k]+1] for k in range(self._Nb)]

    def flow_external(self, uinf, soln, X, Y):
        '''
        Compute net flow U,V on mesh X,Y due to onset flow and all bodies.
        '''
        # Initialize to onset flow
        U,V = uinf[0]*np.ones(X.shape), uinf[1]*np.ones(Y.shape)
        sigma, gamma = soln[0], soln[1]
        for k,body in enumerate(self._bodies):
            # Contribution from source distribution along bodies
            u, v = velocity_source_body(body, sigma[k], X, Y)
            U += u
            V += v
            # Contribution from vorticity distribution along bodies
            u, v = velocity_vortex_body(body, gamma[k]*np.ones(body.nedge),
                X, Y)
            U += u
            V += v
        return U,V

    def pressure_self(self, uinf, soln):
        '''
        Compute pressure coefficient (via Bernoulli) around each body.
        '''
        qt = self.flow_self(uinf, soln)
        return [1. - (qt[k]/nla.norm(uinf))**2 for k in range(self._Nb)]

    def pressure_from_flow(self, uinf, soln, U, V):
        '''
        Compute field of pressure coefficients (via Bernoulli) from flow data.
        '''
        return 1. - (U**2 + V**2)/nla.norm(uinf)**2

    def pressure(self, uinf, soln, X, Y):
        '''
        Compute pressure field on a mesh X,Y from solution.
        '''
        U,V = self.flow_external(uinf, soln, X, Y)
        return self.pressure_from_flow(uinf, soln, U, V)
