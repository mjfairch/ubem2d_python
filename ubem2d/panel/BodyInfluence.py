from ubem2d.panel.PanelInfluence import sf_source_panel
from ubem2d.panel.PanelInfluence import sf_vortex_panel
from ubem2d.panel.PanelInfluence import velocity_source_panel
from ubem2d.panel.PanelInfluence import velocity_vortex_panel
from ubem2d.panel.PanelInfluence import source_influence_matrices
from ubem2d.panel.PanelInfluence import vortex_influence_matrices

__all__ = ['sf_source_body', 'sf_vortex_body', 'velocity_source_body',
    'velocity_vortex_body', 'source_influence_matrices_body', 
    'vortex_influence_matrices_body']

def sf_source_body(body,s,X,Y,m=5):
    '''
    Return the stream function at X,Y due to source sheets of strength s along
    the given body.  Use m subintervals in each Riemann sum.
    '''
    return sf_source_panel(body.x[:-1],body.y[:-1],body.tx,body.ty,body.edge,
        s,X,Y,m)

def sf_vortex_body(body,s,X,Y,m=5):
    '''
    Return the stream function at X,Y due to vortex sheets of strength s along
    the given body.  Use m subintervals in each Riemann sum.
    '''
    return sf_vortex_panel(body.x[:-1],body.y[:-1],body.tx,body.ty,body.edge,
        s,X,Y,m)

def velocity_source_body(body,s,X,Y):
    '''
    Return the velocity at X,Y due to source sheets of strength s along the 
    given body.
    '''
    return velocity_source_panel(body.x[:-1],body.y[:-1],body.tx,body.ty,
        body.edge,s,X,Y)

def velocity_vortex_body(body,s,X,Y):
    '''
    Return the velocity at X,Y due to vortex sheets of strength s along the 
    given body.
    '''
    return velocity_vortex_panel(body.x[:-1],body.y[:-1],body.tx,body.ty,
        body.edge,s,X,Y)

def source_influence_matrices_body(body):
    '''
    Return the self-influence matrices for unit source sheets along the
    given body.
    '''
    return source_influence_matrices(body.x[:-1], body.y[:-1], body.tx,
        body.ty, body.nx, body.ny, body.edge)

def vortex_influence_matrices_body(body):
    '''
    Return the self-influence matrices for unit vortex sheets along the
    given body.
    '''
    return vortex_influence_matrices(body.x[:-1], body.y[:-1], body.tx,
        body.ty, body.nx, body.ny, body.edge)
