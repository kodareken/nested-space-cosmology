"""Directed group22 tail integration with the unchanged order16/ad4 recipe."""
import ast
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import inspect

import mpmath as mp
import sympy as sp

from .nsc_incoming_source_tail import _mp_bloch, leading_reference_matching_identity
from .nsc_incoming_source_quadrature_bound import (
    interval_bloch, gauss_rule, certify_gauss_brackets, pack, unpack, _conj_coeff,
)
from .nsc_incoming_projector_energy_bound import incoming_riccati_interval_coefficients
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range, _up_float


def _root(value): return value**mp.iv.mpf('0.5')


def _gap_root(E, C, B):
    return E/_root(C)*_root(1+C*B/(E*E))


def integer_power(value, exponent):
    """Exact integer algebra; reciprocal BEFORE powers avoids interval wrapping."""
    if not isinstance(exponent,int) or abs(exponent)>32:
        raise ValueError('bounded exact integer power required')
    if exponent < 0:
        value = 1/value; exponent = -exponent
    result = 1
    while exponent:
        if exponent & 1: result = result*value
        exponent >>= 1
        if exponent: value = value*value
    return result


def _literal(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return Fraction(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub): return -_literal(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div): return _literal(node.left)/_literal(node.right)
    raise ValueError('literal integer/half-integer exponent required')


def adapt_bloch_ast(source):
    """Fail-closed guard of generated gap definitions and exact replacements."""
    tree = ast.parse(source); fn = tree.body[0]
    if not isinstance(fn, ast.FunctionDef) or [a.arg for a in fn.args.args] != ['q','mass','angular','frequency']:
        raise ValueError('generated Bloch signature changed')
    assignments = {n.targets[0].id: n.value for n in fn.body if isinstance(n, ast.Assign) and isinstance(n.targets[0], ast.Name)}
    expected = {'x0':'angular**2','x1':'mass**2','x2':'sin(q)','x3':'x2**2',
                'x4':'x3**(-1.0)','x5':'x1*x4','x6':'frequency**2',
                'x7':'2*q','x8':'sin(x7)','x9':'-3*q-x3+(3/2)*x8+3*pi',
                'x10':'x6/x9','x11':'x0+x10+x5'}
    for name, expression in expected.items():
        if name not in assignments or ast.dump(assignments[name]) != ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('generated gap assignment changed: '+name)
    dependent = {'frequency'}
    for name,value in assignments.items():
        if any(isinstance(n,ast.Name) and n.id in dependent for n in ast.walk(value)):
            dependent.add(name)
    for node in ast.walk(tree):
        if isinstance(node,ast.Call) and any(isinstance(n,ast.Name) and n.id in dependent
                                           for arg in node.args for n in ast.walk(arg)):
            allowed = (isinstance(node.func,ast.Name) and (node.func.id == 'array' or
                       node.func.id == 'sqrt' and len(node.args) == 1 and
                       isinstance(node.args[0],ast.Name) and node.args[0].id == 'x11'))
            if not allowed: raise ValueError('unowned energy-dependent analytic function')
        if isinstance(node,ast.BinOp) and isinstance(node.op,ast.Pow) and any(
                isinstance(n,ast.Name) and n.id in dependent for n in ast.walk(node.left)):
            power = _literal(node.right)
            if power.denominator != 1 and not (isinstance(node.left,ast.Name) and node.left.id == 'x11'):
                raise ValueError('unowned energy-dependent fractional power')
    counts = {'square_roots': 0, 'odd_gap_root_powers': []}
    def gap(): return ast.parse('tail_gap_root(frequency,x9,x0+x5)',mode='eval').body
    class Replace(ast.NodeTransformer):
        def visit_Call(self, node):
            node = self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'sqrt' and len(node.args) == 1 and isinstance(node.args[0],ast.Name) and node.args[0].id == 'x11':
                counts['square_roots'] += 1
                return ast.copy_location(gap(), node)
            return node
        def visit_BinOp(self, node):
            node = self.generic_visit(node)
            if isinstance(node.op,ast.Pow) and isinstance(node.left,ast.Name) and node.left.id == 'x11':
                power = _literal(node.right)
                if power.denominator == 2:
                    counts['odd_gap_root_powers'].append(power.numerator)
                    return ast.copy_location(ast.BinOp(gap(),ast.Pow(),ast.Constant(power.numerator)),node)
                if power.denominator != 1: raise ValueError('unowned fractional gap power')
            return node
    tree = Replace().visit(tree)
    if counts != {'square_roots':1,'odd_gap_root_powers':[-3,-5,-7,9]}:
        raise ValueError('generated gap replacement whitelist changed')
    class IntegerAlgebra(ast.NodeTransformer):
        def visit_BinOp(self,node):
            node = self.generic_visit(node)
            if isinstance(node.op,ast.Pow) and any(isinstance(n,ast.Name) and n.id in dependent for n in ast.walk(node.left)):
                exponent = _literal(node.right)
                if exponent.denominator != 1: raise ValueError('gap adapter left a fractional energy power')
                return ast.copy_location(ast.Call(ast.Name('tail_integer_power',ast.Load()),[node.left,ast.Constant(int(exponent))],[]),node)
            if (isinstance(node.op,ast.Div) and isinstance(node.right,ast.Call)
                    and isinstance(node.right.func,ast.Name) and node.right.func.id == 'tail_integer_power'):
                base,exponent = node.right.args
                reciprocal = ast.Call(ast.Name('tail_integer_power',ast.Load()),[base,ast.Constant(-int(_literal(exponent)))],[])
                return ast.copy_location(ast.BinOp(node.left,ast.Mult(),reciprocal),node)
            return node
    tree = IntegerAlgebra().visit(tree); ast.fix_missing_locations(tree)
    return tree, counts


@lru_cache(maxsize=1)
def analytic_bloch():
    source = inspect.getsource(_mp_bloch()); tree, counts = adapt_bloch_ast(source)
    environment = dict(interval_bloch().__globals__)
    environment.update(tail_gap_root=_gap_root,tail_integer_power=integer_power)
    exec(compile(tree,'<owned-tail-gap-branch>','exec'),environment)
    function = environment[tree.body[0].name]
    function.adapter_metadata = {'generated_source_sha256':sha256(source.encode()).hexdigest(),
                                 'replacement_whitelist':counts,
                                 'integer_power_policy':'integer multiplication; reciprocal base before powers; divide by power rewritten as multiply by inverse power',
                                 'branch':'E/sqrt(C)*sqrt(1+C*B/E^2); positive real agreement'}
    return function


@lru_cache(maxsize=1)
def adiabatic_grade_certificate():
    """Propagate x-orders through the already-owned b0..b4 recurrence.

    q differentiation preserves powers. Its only special leading case is
    b0_z=1+O(x^2), whose metric-independent constant differentiates to zero.
    h_transverse=O(1), h_z=O(x^-1), inverse_energy_squared=O(x^2).
    """
    grades = [(1,0)]; rows = []
    for n in range(1,5):
        transverse,longitudinal = (1,2) if n == 1 else grades[-1]
        cross = (min(longitudinal+2,transverse+1),transverse+2)
        normal = min([min(grades[i][0]+grades[n-i][0],grades[i][1]+grades[n-i][1])
                      for i in range(1,n)] or [1000])
        actual = (min(cross[0],normal+1),min(cross[1],normal))
        if actual != (n+1,n+2): raise ArithmeticError('owned adiabatic recurrence power bound failed')
        grades.append(actual)
        reduced = min(actual[0]-2,actual[1]-3)
        rows.append({'order':n,'transverse_power':actual[0],'longitudinal_power':actual[1],
                     'source_over_x_squared_minimum_power':reduced,
                     'grade_residuals':[actual[0]-(n+1),actual[1]-(n+2)]})
    return {'owner':'nsc_compact_ctp_neck._adiabatic_bloch_function',
            'recurrence':'b_n=-h cross d(b_(n-1))/(2 energy^2)-sum_(i=1..n-1)(b_i dot b_(n-i))*b0/2',
            'leading_b0':'transverse O(x); longitudinal 1+O(x^2), leading1 independent of q',
            'rows':rows,'ad2_ad4_have_no_remaining_poles':all(r['source_over_x_squared_minimum_power']>=1 for r in rows[1:])}


@lru_cache(maxsize=1)
def endpoint_certificate():
    """Exact rearrangement plus the reused c1/c2/ad1 removability owner."""
    root, y, x, U, R = sp.symbols('root y x U R')
    base = U*x/(1+root); difference = R+U*x*y/(2*(1+root)**2)
    residual = sp.factor((base+difference-U*x/2-R).subs(y,root*root-1))
    norm = y/(1+root)**2
    longitudinal = sp.factor(((1-norm)/(1+norm)-1/root).subs(y,root*root-1))
    transverse = sp.factor((2/((1+root)*(1+norm))-1/root).subs(y,root*root-1))
    inherited = leading_reference_matching_identity()
    if [residual,longitudinal,transverse] != [0,0,0] or inherited['vacuum_minus_ad1'] != ['0']*3:
        raise ArithmeticError('exact endpoint cancellation failed')
    return {'stable_rearrangement_residuals':[str(v) for v in (residual,longitudinal,transverse)],
            'reused_c1_c2_ad1': inherited,
            'higher_adiabatic_grade_certificate':adiabatic_grade_certificate(),
            'meromorphic_to_analytic':'only finite integer x poles after root factoring; owned source=x^2*g removes them',
            'endpoint_g0':0, 'numerical_chopping':False}


class CompactTailKernel:
    def __init__(self, channel, sign):
        if channel['index'] != 22 or sign not in (1,-1): raise ValueError('bounded group22 signed pilot required')
        self.a,self.r,self.mass,angular,self.factor = _parameters(channel)
        self.factor /= 2; self.ell = sign*angular
        positive = incoming_riccati_interval_coefficients(channel,order=16)
        self.c = positive if sign == 1 else [(-1)**(j+1)*_conj_coeff(v) for j,v in enumerate(positive)]
        self.cb = [_conj_coeff(v) for v in self.c]
        self.bloch = analytic_bloch()

    def margins(self, max_x):
        q = 3*mp.iv.pi/4
        C = 3*(mp.iv.pi-q)+mp.iv.mpf('1.5')*mp.iv.sin(2*q)-mp.iv.sin(q)**2
        B = self.mass**2/mp.iv.sin(q)**2+self.ell**2
        gap = C*B*max_x**2
        normalization = self.a**2*sum((abs(c)*(max_x/2)**(j+1) for j,c in enumerate(self.c)),mp.iv.mpf(0))**2
        if _hi(gap) >= 1 or _hi(normalization) >= 1: raise ArithmeticError('tail disk pole/branch exclusion failed')
        return {'gap_product_upper':pack(gap),'normalization_product_upper':pack(normalization)}

    def __call__(self, x):
        a,r,m,ell = self.a,self.r,self.mass,self.ell
        U,Ub = ell/r+mp.iv.j*m,ell/r-mp.iv.j*m
        y = a*a*(m*m+(ell/r)**2)*x*x; root = _root(1+y)
        base,barbase = U*x/(1+root),Ub*x/(1+root)
        remainder = sum((self.c[j]*(x/2)**(j+1) for j in range(1,16)),mp.iv.mpc(0))
        barremainder = sum((self.cb[j]*(x/2)**(j+1) for j in range(1,16)),mp.iv.mpc(0))
        dS = remainder+U*x*y/(2*(1+root)**2)
        dSb = barremainder+Ub*x*y/(2*(1+root)**2)
        D = 1+a*a*base*barbase
        du = a*a*(base*dSb+barbase*dS+dS*dSb)
        denominator = D*(D+du)
        reS,imS = (base+barbase)/2,(base-barbase)/(2*mp.iv.j)
        reD,imD = (dS+dSb)/2,(dS-dSb)/(2*mp.iv.j)
        delta = [2*a*(D*imD-imS*du)/denominator,
                 -2*a*(D*reD-reS*du)/denominator,-2*du/denominator]
        terms = self.bloch(3*mp.iv.pi/4,m,ell,1/x)
        for i in range(3): delta[i] -= sum((terms[j][i,0] for j in range(1,5)),mp.iv.mpf(0))
        parallel = -delta[2]/(a*x)
        return [(-m*delta[0]+ell/r*delta[1]+parallel)/(x*x),
                parallel/(x*x),mp.iv.mpc(0),ell*delta[1]/(2*r*x*x)]


def prepare_integral(channel, *, precision=80):
    endpoint = endpoint_certificate()
    with _precision(precision):
        rule = gauss_rule(48,precision); nodes = certify_gauss_brackets(rule['roots'],48)
        h = mp.iv.mpf(1)/320; radius = 2*h; max_x = 3*h
        signs = []
        for sign in (1,-1):
            kernel = CompactTailKernel(channel,sign); margins = kernel.margins(max_x)
            maxima = [mp.mpf(0)]*4
            for arc in range(32):
                theta = 2*mp.iv.pi*_range(mp.mpf(arc)/32,mp.mpf(arc+1)/32)
                x = mp.iv.mpc(h+radius*mp.iv.cos(theta),radius*mp.iv.sin(theta))
                for j,value in enumerate(kernel(x)):
                    maxima[j] = max(maxima[j],_hi(abs(value)))
            if not all(mp.isfinite(v) for v in maxima): raise ArithmeticError('nonfinite analytic tail circle enclosure')
            samples = []
            for node,_ in nodes:
                values = kernel(h+h*node)
                if any(not _lo(v.imag) <= 0 <= _hi(v.imag) for v in values): raise ArithmeticError('real tail sample lost reality')
                samples.append([pack(v.real) for v in values])
            signs.append({'sign':sign,'physical_factor':pack(kernel.factor),'margins':margins,
                          'circle_absolute_upper':[pack(mp.iv.mpf(v)) for v in maxima], 'samples':samples})
            print(f'group22 order16 tail: signed integral {sign} prepared',flush=True)
        return {'group':22,'lower_energy':160.,'riccati_order':16,'reference_order':4,
                'precision':precision,'points':48,'circle_arcs':32,'halfwidth':pack(h),'radius':pack(radius),
                'rule':rule,'signs':signs,'endpoint_certificate':endpoint,
                'analytic_adapter':analytic_bloch().adapter_metadata}


def replay_integral(payload):
    if (payload['group'],payload['lower_energy'],payload['riccati_order'],payload['reference_order'],payload['points']) != (22,160.,16,4,48):
        raise ValueError('unchanged group22 order16/ad4 tail and fixed48 rule required')
    if payload['endpoint_certificate'] != endpoint_certificate(): raise ValueError('endpoint identity changed')
    with _precision(payload['precision']):
        nodes = certify_gauss_brackets(payload['rule']['roots'],48)
        if [pack(w) for _,w in nodes] != payload['rule']['weights']: raise ValueError('certified weights changed')
        h,radius = unpack(payload['halfwidth']),unpack(payload['radius'])
        if pack(h) != pack(mp.iv.mpf(1)/320) or pack(radius) != pack(2*h): raise ValueError('fixed compactified cell changed')
        total = [mp.iv.mpf(0) for _ in range(4)]; errors = [mp.iv.mpf(0) for _ in range(4)]
        if [r['sign'] for r in payload['signs']] != [1,-1]: raise ValueError('actual angular signs required')
        for row in payload['signs']:
            if len(row['samples']) != 48 or any(_hi(unpack(v)) >= 1 for v in row['margins'].values()):
                raise ValueError('sample coverage or pole margin failed')
            factor = unpack(row['physical_factor'])
            for j in range(4):
                quadrature = h*sum((w*unpack(sample[j]) for (_,w),sample in zip(nodes,row['samples'])),mp.iv.mpf(0))
                M = unpack(row['circle_absolute_upper'][j]); error = 4*h*M*(h/radius)**96/(1-h/radius)
                total[j] += factor*quadrature; errors[j] += factor*error
        enclosed = [v+_range(-_hi(e),_hi(e)) for v,e in zip(total,errors)]
        return {'integral_interval':[pack(v) for v in enclosed],
                'quadrature_interval':[pack(v) for v in total],
                'analytic_quadrature_error_upper':[pack(v) for v in errors]}


def directed_thermal_tail(channel,config):
    """Same elementary tail formula; outward intervals preserve nonzero values."""
    with _precision(80):
        a,r,m,ell,factor = _parameters(channel)
        kappa,omega = mp.iv.mpf(config['surface_gravity']),mp.iv.mpf(config['omega'])
        m0,m1 = mp.iv.mpf(0),mp.iv.mpf(0)
        for coefficient,alpha in ((2,mp.iv.pi/kappa),(1,2*mp.iv.pi/(omega*kappa))):
            exponential = coefficient*mp.iv.exp(-alpha*160)
            m0 += exponential/alpha; m1 += exponential*(160/alpha+1/alpha**2)
        M = mp.iv.sqrt(m*m+(ell/r)**2)
        values = [factor*(m1/a+M*m0),factor*m1/a,factor*m1/a,factor*ell*m0/(2*r)]
        if any(_lo(v) <= 0 for v in values): raise ArithmeticError('nonzero group22 thermal tail lost')
        return {'intervals':[pack(v) for v in values],
                'upper_display':[mp.nstr(_hi(v),25) for v in values],
                'display_is_not_proof':'lossless interval endpoints own outward rounding',
                'source_law':'2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
                'exact_zero_components':[]}
