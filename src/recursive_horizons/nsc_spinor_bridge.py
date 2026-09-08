"""Exact Lorentz/chiral/charge checks for a candidate two-sheet mass bridge.

The scalar sheet bilinear is a kinematic benchmark, not a derived throat map.
Sheet, angular-channel, radial and physical chiral labels remain distinct.
"""
import sympy as sp
import numpy as np
from scipy.linalg import expm, eigh


def pauli():
    return (sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))


def weyl_matrices():
    s=pauli();i=sp.eye(2);z=sp.zeros(2)
    g0=z.row_join(i).col_join(i.row_join(z))
    gi=[z.row_join(a).col_join((-a).row_join(z)) for a in s]
    g5=sp.diag(-1,-1,1,1)
    return {'gamma':[g0,*gi],'gamma5':g5,'left':(sp.eye(4)-g5)/2,
            'right':(sp.eye(4)+g5)/2,'alpha':[g0*a for a in gi],
            'parity':g0,'charge_conjugation_matrix':sp.I*gi[1]}


def comm(a,b):return a*b-b*a
def anti(a,b):return a*b+b*a


def exact_checks():
    w=weyl_matrices();g=w['gamma'];g5=w['gamma5'];i4=sp.eye(4);s=pauli()
    px,py,pz,m,phi,e,z,k,r,v=sp.symbols('p_x p_y p_z m Phi charge z kappa radius vertex',real=True)
    p=[px,py,pz];p2=sum(a*a for a in p)
    kinetic=sum((a*b for a,b in zip(p,w['alpha'])),sp.zeros(4))
    h=kinetic+m*g[0];bc=w['charge_conjugation_matrix']
    checks={}
    def zero(name,value):
        value=sp.simplify(value)
        if value != sp.zeros(*value.shape):raise AssertionError((name,value))
        checks[name]='0'
    for a in range(4):
        for b in range(4):
            metric=(1 if a==0 else -1) if a==b else 0
            zero(f'Clifford_{a}_{b}',anti(g[a],g[b])-2*metric*i4)
    zero('projectors_resolve_identity',w['left']+w['right']-i4)
    zero('projectors_are_orthogonal',w['left']*w['right'])
    zero('free_H_squared',h*h-(p2+m*m)*i4)
    zero('massless_chirality_conservation',comm(kinetic,g5))
    zero('mass_chiral_commutator',comm(h,g5)-2*m*g[0]*g5)
    zero('parity_flips_chirality',g[0]*g5*g[0]+g5)
    zero('parity_reverses_momentum',g[0]*h*g[0]-(-kinetic+m*g[0]))
    zero('C_cnumber_flips_chirality',bc*g5.conjugate()*bc.H+g5)
    zero('C_cnumber_reverses_frequency',bc*h.conjugate()*bc.H+(-kinetic+m*g[0]))
    zero('C_maps_L_to_R_conjugate',bc*w['left'].conjugate()-w['right']*bc)

    # Tensor order sheet (tau) then full Weyl spinor. The mass in Hamiltonian
    # includes beta: a scalar bar(psi_p) psi_c becomes psi_p† beta psi_c.
    i2=sp.eye(2);kp=sp.kronecker_product
    full=kp(i2,kinetic)+phi*kp(s[0],g[0])
    gamma5=kp(i2,g5);sheet=kp(s[2],i4)
    zero('scalar_sheet_mass_gap',full*full-(p2+phi*phi)*sp.eye(8))
    zero('sheet_grading_commutes_with_physical_chirality',comm(sheet,gamma5))
    zero('scalar_bridge_chiral_mixing',comm(full,gamma5)-2*phi*kp(s[0],g[0]*g5))
    # Exact Schur self energy phi² beta(z-alpha.p)^-1 beta.
    inverse_child=(z*i4+kinetic)/(z*z-p2)
    zero('scalar_bridge_shared_self_energy',phi**2*g[0]*inverse_child*g[0]-phi**2*(z*i4-kinetic)/(z*z-p2))
    pseudo=sp.I*g[0]*g5
    zero('pseudoscalar_H_vertex_is_Hermitian',pseudo-pseudo.H)
    pseudo_full=kp(i2,kinetic)+phi*kp(s[0],pseudo)
    zero('pseudoscalar_has_same_gap',pseudo_full*pseudo_full-(p2+phi*phi)*sp.eye(8))
    zero('pseudoscalar_has_same_scalar_self_energy',phi**2*pseudo*inverse_child*pseudo-phi**2*(z*i4-kinetic)/(z*z-p2))
    naive=kp(i2,kinetic)+phi*kp(s[0],i4)
    zero('identity_sheet_link_is_not_mass_gap',naive*naive-(p2+phi*phi)*sp.eye(8)-2*phi*kp(s[0],kinetic))
    zero('same_charge_sheet_U1',comm(e*sp.eye(8),full))
    zero('opposite_charge_sheet_breaking',comm(e*sheet,full)-2*sp.I*e*phi*kp(s[1],g[0]))

    # A constructive restriction: parent-left plus child-right is invariant
    # for THIS scalar coupling. Its equally invariant complement is not
    # excluded by the operator; a physical selection has not been derived.
    projector=(sp.eye(8)-sheet*gamma5)/2
    embedding=sp.eye(8)[:,[0,1,6,7]]
    zero('restricted_sector_embedding',embedding*embedding.H-projector)
    zero('restricted_sector_preserved',comm(projector,full))
    zero('restricted_sector_is_one_Dirac_H',embedding.H*full*embedding-(kinetic+phi*g[0]))
    zero('restricted_joint_parent_R_vanishes',kp((i2+s[2])/2,w['right'])*projector)
    zero('restricted_joint_child_L_vanishes',kp((i2-s[2])/2,w['left'])*projector)
    paired_parity=kp(s[0],g[0])
    zero('combined_sheet_parity_preserves_restriction',comm(paired_parity,projector))
    zero('combined_sheet_parity_restricts_to_Dirac_parity',embedding.H*paired_parity*embedding-g[0])
    swap=kp(s[0],i4)
    zero('sheet_swap_exchanges_complementary_restrictions',swap*projector*swap-(sp.eye(8)-projector))
    paired_c=kp(s[0],bc)
    zero('combined_sheet_C_preserves_restriction',paired_c*projector.conjugate()*paired_c.H-projector)
    zero('combined_sheet_C_restricts_to_Dirac_C',embedding.H*paired_c*embedding.conjugate()-bc)

    # Angular eta, radial rho; these eta matrices are NOT geometric sheets.
    radial=kp(i2,pz*s[1])+kp(s[2],k/r*s[0])
    physical_g5=-kp(s[0],s[1]);beta=kp(i2,s[2])
    pulse=kp(s[2],v*s[0])
    zero('paired_angular_chirality_squared',physical_g5*physical_g5-i4)
    zero('paired_massless_radial_chirality',comm(radial,physical_g5))
    zero('radius_pulse_preserves_massless_chirality',comm(pulse,physical_g5))
    zero('physical_radial_mass_breaks_chirality',comm(radial+m*beta,physical_g5)-2*m*beta*physical_g5)
    zero('single_radial_grading_reflects_energy',anti(s[2],pz*s[1]+k/r*s[0]))
    normal_reflection=kp(s[0],s[2])
    zero('symmetric_throat_normal_reflection',normal_reflection*(radial.subs(pz,-pz)+m*beta)*normal_reflection-(radial+m*beta))
    zero('throat_reflection_flips_chirality',normal_reflection*physical_g5*normal_reflection+physical_g5)
    ranks={f'sheet_{a}_chirality_{b}':int((kp((i2+a*s[2])/2,(i4+b*g5)/2)).rank())
           for a in (-1,1) for b in (-1,1)}
    return {'identities':checks,'joint_projector_ranks':ranks,
            'candidate_restriction':{'projector':'(I8-tau3 tensor gamma5)/2',
                'selected_components':'parent-left and child-right','rank':int(projector.rank()),
                'complement_rank':int((sp.eye(8)-projector).rank()),
                'embedding_columns':[0,1,6,7],
                'dynamically_invariant_for_scalar_bridge':True,
                'selected_by_actual_throat_or_action':False},
            'matrices':{'gamma5_Weyl':str(g5),'parity':str(g[0]),'C_cnumber_B':str(bc),
                        'paired_angular_gamma5':str(physical_g5)},
            'unresolved_bridge':'derive the actual Clifford structure and any justified local scalar or pseudoscalar limit from the curved boundary response'}


def chiral_evolution(momentum=1.2,mass=.7,helicity=1,times=(0.,.5,1.,2.)):
    """Fixed-helicity benchmark; a pure chiral vector contains both energies."""
    if helicity not in (-1,1):raise ValueError('helicity sign must be +/-1')
    h=np.array([[-helicity*momentum,mass],[mass,helicity*momentum]])
    e,u=eigh(h);energy=float(np.hypot(momentum,mass))
    left=np.array([1.,0.]);positive=u[:,1]
    controls=[]
    for time in times:
        state=expm(-1j*h*time)@left
        evolved=expm(-1j*h*time)@positive
        analytic=(mass/energy)**2*np.sin(energy*time)**2 if energy else 0.
        controls.append({'time':time,'pure_left_to_right':float(abs(state[1])**2),
                         'analytic_transition':float(analytic),
                         'positive_energy_right_probability':float(abs(evolved[1])**2)})
    return {'momentum':momentum,'mass':mass,'helicity_sign':helicity,'energy':energy,
            'pure_left_negative_energy_weight':float(abs(np.vdot(u[:,0],left))**2),
            'positive_energy_chirality':float(abs(positive[1])**2-abs(positive[0])**2),
            'expected_positive_energy_chirality':float(helicity*momentum/energy) if energy else 0.,
            'controls':controls}
