"""Whole supplied-family derivatives, explicit data, and coherent contraction."""
import numpy as np
import pytest
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import expm_multiply

from recursive_horizons.nsc_prepared_history_jets import (
    evolve_prepared_field_jets, restriction_covariance_tangent,
)
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator


class SuppliedControl:
    """Synthetic PG coefficient family; neither a parent nor a KS solution."""
    def __init__(self, owner, parameters, ndir=2):
        self.owner, self.parameters, self.ndir = owner, np.asarray(parameters), ndir

    def shapes(self, time, x):
        bump = np.sin(np.pi*(x-x[0])/(x[-1]-x[0]))**2
        bump[[0, -1]] = 0.
        return (1+.4*time)*bump, (1-.3*time)*bump*(.8+.1*x)

    def values(self, time, x):
        h, j = self.shapes(time, x);p, q = self.parameters
        result = self.owner.reference.copy()
        result[0] *= np.exp(p*h)
        result[1] += q*j
        result[3] *= np.exp(.4*p*h+.6*q*j)
        return result

    def log_directions(self, time, x, metric):
        h, j = self.shapes(time, x)
        result = np.zeros((2, 4, len(x)))
        result[0, 0] = h;result[0, 3] = .4*h
        result[1, 1] = j;result[1, 3] = .6*j
        return result[:self.ndir]


@pytest.fixture
def control():
    x = np.linspace(-1.3, 1.3, 21)
    owner = FourthOrderModePropagator(x, .7, 1.1)
    row = np.arange(2*len(x))[:, None];column = np.arange(3)[None, :]
    base = .2*np.exp(.08j*row*(1+column))
    dphi = np.array([.13*np.cos(.12*row+column)+.04j,
                    -.09j*np.sin(.07*row-column)+.05])
    source = np.zeros_like(base)
    source[[len(x)-1, 2*len(x)-1]] = np.array([[.12+.03j, .07j, -.04], [.02j, .11, .09-.02j]])
    dsource = np.array([(.6+.2j)*source, (-.4+.5j)*source])
    parameters = np.array([.012, -.009]);times = np.linspace(0., .06, 7)
    def initial(parameters):return base+np.einsum('d,drs->rs', parameters, dphi)
    def incoming(parameters, ndir=2):
        def supplied(time):
            forcing=(1+.7*time)*source+np.einsum('d,drs->rs', parameters, dsource)*(1-.4*time)
            return forcing, dsource[:ndir]*(1-.4*time)
        return supplied
    return owner, parameters, times, initial, dphi, incoming


def independent_primal(owner, times, provider, initial, incoming):
    """No tangent construction; directly integrate the whole perturbed family."""
    rows, sources = initial.shape
    y = np.vstack((initial, np.eye(sources)))
    for left, right in zip(times[:-1], times[1:]):
        t = (left+right)/2
        L, _ = owner.generator(provider.values(t, owner.x))
        f, _ = incoming(t)
        operator = bmat([[L, csr_matrix(f)], [None, csr_matrix((sources, sources))]], format='csr')
        dt = right-left
        y = expm_multiply(dt*operator, y, traceA=dt*operator.diagonal().sum())
    return y[:rows]


def test_full_family_tangent_contains_nonzero_preparation_and_incident_terms(control):
    owner, p, times, initial, dphi, incoming = control
    provider = SuppliedControl(owner, p)
    result = evolve_prepared_field_jets(owner, times, provider, initial(p), dphi, incoming(p))
    primal = independent_primal(owner, times, provider, initial(p), incoming(p))
    assert np.max(abs(primal-result['evolved_field'])) < 2e-14
    residuals=[]
    for direction in np.eye(2):
        step=2e-5
        plus,minus=[independent_primal(owner,times,SuppliedControl(owner,p+sign*step*direction),
                    initial(p+sign*step*direction),incoming(p+sign*step*direction)) for sign in (1,-1)]
        index=np.argmax(direction)
        residuals.append(np.max(abs((plus-minus)/(2*step)-result['tangent_fields'][index])))
    assert max(residuals) < 3e-10
    without_initial=evolve_prepared_field_jets(owner,times,provider,initial(p),np.zeros_like(dphi),incoming(p))
    def without_forcing_tangent(t):
        f,df=incoming(p)(t);return f,np.zeros_like(df)
    without_forcing=evolve_prepared_field_jets(owner,times,provider,initial(p),dphi,without_forcing_tangent)
    assert np.max(abs(result['tangent_fields']-without_initial['tangent_fields'])) > .01
    assert np.max(abs(result['tangent_fields']-without_forcing['tangent_fields'])) > 1e-4
    assert result['norm_flux_residual'] < 3e-11
    assert result['metric_flux_tangent_residual'] < 3e-11
    assert result['auxiliary_identity_residual'] < 3e-14
    assert result['physical_preparation_status']==result['physical_history_status']=='OPEN'
    assert result['continuum_error_bound'] is None
    assert result['sampled_causal_buffer'] is None
    assert result['full_branch_unitarity_claimed'] is False
    print('prepared-field residuals', list(map(float,residuals)), 'flux',result['norm_flux_residual'],
          'dflux',result['metric_flux_tangent_residual'],'auxiliary',result['auxiliary_identity_residual'])


def test_zero_directions_remain_explicit_and_match_primal(control):
    owner,p,times,initial,dphi,incoming=control
    provider=SuppliedControl(owner,p,ndir=0)
    result=evolve_prepared_field_jets(owner,times,provider,initial(p),dphi[:0],incoming(p,ndir=0))
    assert result['tangent_fields'].shape==(0,*initial(p).shape)
    assert np.max(abs(result['evolved_field']-independent_primal(owner,times,provider,initial(p),incoming(p)))) < 2e-14


def test_coherent_restriction_keeps_all_three_terms_and_direction_axes():
    F=np.array([[1,.2j,.3],[.1,.7,-.15j]],complex)
    dF=np.array([[[.1j,.02,.03],[.04,.06j,.02]],[[.03,.02j,-.07],[.1j,-.01,.04]]])
    C=np.array([[.6,.08+.04j,.03j],[.08-.04j,.4,.06],[-.03j,.06,.5]])
    dC=np.array([[[.04,.02j,.03],[-.02j,-.01,.01j],[.03,-.01j,.02]],
                 [[-.03,.04,.01j],[.04,.02,.02],[-.01j,.02,.01]]])
    result=restriction_covariance_tangent(F,dF,C,dC)
    errors=[]
    for j in range(2):
        h=2e-5
        values=[(F+sign*h*dF[j])@(C+sign*h*dC[j])@(F+sign*h*dF[j]).conj().T for sign in (1,-1)]
        errors.append(np.max(abs((values[0]-values[1])/(2*h)-result[j])))
        assert np.max(abs(result[j]-F@dC[j]@F.conj().T)) > .01
        assert np.max(abs(result[j]-restriction_covariance_tangent(F,dF[j],C,np.zeros_like(C)))) > .01
    assert max(errors) < 3e-11
    assert np.max(abs(result-result.swapaxes(-1,-2).conj())) < 3e-15
    # Preserve coherent terms involving the third source column; never silently truncate it.
    truncated=restriction_covariance_tangent(F[:,:2],dF[:,:,:2],C[:2,:2],dC[:,:2,:2])
    assert np.max(abs(result-truncated)) > .005
    print('coherent restriction residuals',list(map(float,errors)))


@pytest.mark.parametrize('times',[None,[0.],[0.,0.],[.1,0.],[0.,np.nan],[0.,np.inf],[0.,.1j]])
def test_invalid_time_coordinates_rejected(control,times):
    owner,p,_,initial,dphi,incoming=control
    with pytest.raises(ValueError):
        evolve_prepared_field_jets(owner,times,SuppliedControl(owner,p),initial(p),dphi,incoming(p))


@pytest.mark.parametrize('missing',['initial','incoming','incoming_tangent','metric_tangent','direction_count'])
def test_missing_preparation_inputs_never_become_implicit_zero(control,missing):
    owner,p,times,initial,dphi,incoming=control
    provider=SuppliedControl(owner,p);callback=incoming(p);tangent=dphi
    if missing=='initial':tangent=None
    if missing=='incoming':callback=None
    if missing=='incoming_tangent':callback=lambda t:(np.zeros_like(initial(p)),None)
    if missing=='metric_tangent':provider.log_directions=lambda *args:None
    if missing=='direction_count':provider.ndir=1
    with pytest.raises(ValueError):
        evolve_prepared_field_jets(owner,times,provider,initial(p),tangent,callback)


def test_covariance_missing_or_mismatched_terms_rejected():
    with pytest.raises(ValueError,match='explicit'):
        restriction_covariance_tangent(np.eye(2),np.eye(2),np.eye(2),None)
    with pytest.raises(ValueError,match='dimensions'):
        restriction_covariance_tangent(np.ones((2,3)),np.ones((2,3)),np.eye(2),np.eye(2))
    with pytest.raises(ValueError,match='leading'):
        restriction_covariance_tangent(np.ones((2,2,3)),np.ones((3,2,3)),np.eye(3),np.eye(3))


@pytest.mark.parametrize('kind,row',[('forcing',0),('forcing',4),('tangent',0),('tangent',4)])
def test_incident_callback_rejects_bulk_or_left_boundary_source(control,kind,row):
    owner,p,times,initial,dphi,incoming=control
    def wrong(time):
        f,df=incoming(p)(time)
        if kind=='forcing':f[row,0]=1e-20
        else:df[0,row,0]=1e-20
        return f,df
    with pytest.raises(ValueError,match='right SAT inflow rows'):
        evolve_prepared_field_jets(owner,times,SuppliedControl(owner,p),initial(p),dphi,wrong)


def test_requested_time_trace_samples_preserve_initial_final_and_evolution(control):
    owner,p,times,initial,dphi,incoming=control
    provider=SuppliedControl(owner,p)
    rows=np.array([len(owner.x)-2,2*len(owner.x)-2],dtype=np.int64)
    plain=evolve_prepared_field_jets(owner,times,provider,initial(p),dphi,incoming(p))
    sampled=evolve_prepared_field_jets(owner,times,provider,initial(p),dphi,incoming(p),sample_rows=rows)
    assert not {'sample_times','sample_rows','sampled_fields','sampled_tangent_fields'} & plain.keys()
    assert np.array_equal(sampled['evolved_field'],plain['evolved_field'])
    assert np.array_equal(sampled['tangent_fields'],plain['tangent_fields'])
    assert np.array_equal(sampled['sample_times'],times)
    assert np.array_equal(sampled['sample_rows'],rows)
    assert sampled['sampled_fields'].shape==(len(times),2,initial(p).shape[1])
    assert sampled['sampled_tangent_fields'].shape==(2,len(times),2,initial(p).shape[1])
    assert np.array_equal(sampled['sampled_fields'][0],initial(p)[rows])
    assert np.array_equal(sampled['sampled_tangent_fields'][:,0],dphi[:,rows])
    assert np.array_equal(sampled['sampled_fields'][-1],sampled['evolved_field'][rows])
    assert np.array_equal(sampled['sampled_tangent_fields'][:,-1],sampled['tangent_fields'][:,rows])
    # The interior node is a completed prefix evolution, never interpolation.
    prefix=evolve_prepared_field_jets(owner,times[:4],provider,initial(p),dphi,incoming(p))
    assert np.array_equal(sampled['sampled_fields'][3],prefix['evolved_field'][rows])
    assert np.array_equal(sampled['sampled_tangent_fields'][:,3],prefix['tangent_fields'][:,rows])
    zero=evolve_prepared_field_jets(owner,times,SuppliedControl(owner,p,ndir=0),initial(p),dphi[:0],incoming(p,ndir=0),sample_rows=rows)
    assert zero['sampled_tangent_fields'].shape==(0,len(times),2,initial(p).shape[1])


@pytest.mark.parametrize('rows',[[1.0,2.0],[1,1],[-1],[42],[True,False],[],1,[[1,2]]])
def test_invalid_sample_rows_rejected_before_evolution(control,rows):
    owner,p,times,initial,dphi,incoming=control
    with pytest.raises(ValueError,match='sample_rows'):
        evolve_prepared_field_jets(owner,times,SuppliedControl(owner,p),initial(p),dphi,incoming(p),sample_rows=rows)
