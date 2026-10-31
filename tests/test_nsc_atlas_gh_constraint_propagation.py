"""Independent covariant GH-constraint propagation and damping factors.

These are identities on the exact closed de Sitter reference, rather than
acceptance checks for a truncated nonlinear evolution. In particular a late
decaying homogeneous constraint mode does not bound unresolved spatial modes.
"""
import sympy as sp


def test_trace_reversal_removes_only_the_metric_damping_term():
    gamma=sp.symbols('gamma',real=True)
    C=sp.Matrix(sp.symbols('C0:4',real=True))
    n=sp.Matrix([-1,0,0,0]);g=sp.diag(-1,1,1,1);inv=g.inv()
    nC=(n.T*inv*C)[0]
    damping=gamma*(n*C.T+C*n.T-g*nC)/2
    trace=sp.trace(inv*damping)
    assert sp.simplify(trace+gamma*nC)==0
    assert sp.simplify(damping-g*trace/2-gamma*(n*C.T+C*n.T)/2)==sp.zeros(4)
    # Thus the Bianchi divergence, after multiplying by two, carries
    # -gamma*div(n_a C_b+n_b C_a), not half or twice this coefficient.


def test_homogeneous_covector_equation_from_closed_flrw_connections():
    t,x,y,z=sp.symbols('t x y z',real=True)
    coordinates=(t,x,y,z)
    a=sp.Function('a')(t);c=sp.Function('c')(t)
    h=sp.diff(a,t)/a
    g=sp.diag(-1,a*a,a*a*sp.sin(x)**2,a*a*sp.sin(x)**2*sp.sin(y)**2)
    inv=g.inv()
    connection=sp.MutableDenseNDimArray.zeros(4,4,4)
    for upper in range(4):
        for i in range(4):
            for j in range(4):
                connection[upper,i,j]=sp.simplify(sum(inv[upper,k]*(
                    sp.diff(g[k,j],coordinates[i])+sp.diff(g[k,i],coordinates[j])
                    -sp.diff(g[i,j],coordinates[k]))/2 for k in range(4)))

    C=(c,0,0,0);n=(-1,0,0,0)
    DC=sp.MutableDenseNDimArray.zeros(4,4)
    for i in range(4):
        for j in range(4):
            DC[i,j]=sp.diff(C[j],coordinates[i])-sum(connection[k,i,j]*C[k] for k in range(4))
    boxC0=sp.simplify(sum(inv[i,j]*(sp.diff(DC[j,0],coordinates[i])
        -sum(connection[k,i,j]*DC[k,0]+connection[k,i,0]*DC[j,k] for k in range(4)))
        for i in range(4) for j in range(4)))
    assert sp.simplify(boxC0+sp.diff(c,t,2)+3*h*sp.diff(c,t)-3*h*h*c)==0

    tensor=sp.Matrix(4,4,lambda i,j:n[i]*C[j]+n[j]*C[i])
    damping_div0=sp.simplify(sum(inv[i,j]*(sp.diff(tensor[j,0],coordinates[i])
        -sum(connection[k,i,j]*tensor[k,0]+connection[k,i,0]*tensor[j,k] for k in range(4)))
        for i in range(4) for j in range(4)))
    assert sp.simplify(damping_div0-2*sp.diff(c,t)-6*h*c)==0

    Ric00=sp.simplify(sum(sp.diff(connection[k,0,0],coordinates[k])
        -sp.diff(connection[k,k,0],t)
        +sum(connection[k,k,j]*connection[j,0,0]-connection[k,0,j]*connection[j,k,0]
             for j in range(4)) for k in range(4)))
    assert sp.simplify(Ric00+3*sp.diff(a,t,2)/a)==0
    # R_0{}^a C_a = -R_00*c in the (-+++) signature.
    gamma,H=sp.symbols('gamma H',positive=True)
    propagated=sp.expand(-boxC0+Ric00*c+gamma*damping_div0)
    expected=(sp.diff(c,t,2)+(3*h+2*gamma)*sp.diff(c,t)
              +(6*gamma*h-3*(h*h+H*H))*c)
    assert sp.simplify((propagated-expected).subs(sp.diff(a,t,2),H*H*a))==0


def test_current_and_stronger_damping_late_characteristic_roots():
    H=sp.symbols('H',positive=True)
    x=sp.symbols('x',real=True)  # h/H=tanh(Ht), hence 0<=x<1 on this branch.
    rate=sp.symbols('rate',real=True)
    coefficient=lambda gamma:6*gamma*H*x-3*H*H*(x*x+1)
    assert sp.factor(coefficient(H))==-3*H*H*(x-1)**2
    assert sp.factor(rate**2+5*H*rate+coefficient(H).subs(x,1))==rate*(rate+5*H)
    assert sp.factor(rate**2+7*H*rate+coefficient(2*H).subs(x,1))==(rate+H)*(rate+6*H)
    assert sp.simplify(coefficient(2*H).subs(x,2-sp.sqrt(3)))==0
    # The stronger choice is still initially negative; its positive restoring
    # coefficient begins only after tanh(Ht)>2-sqrt(3), not at the turning slice.
    assert coefficient(2*H).subs(x,0)==-3*H*H
    assert coefficient(2*H).subs(x,1)==6*H*H
