"""Operand-derived binary64 guards for immutable source receipt replay.

These bounds compare two executions of finite weighted-array reductions.
They do not certify quadrature, fields, physical remainders or local gates.
"""
from copy import deepcopy
import numpy as np


def real_values(value):
    """Reject coercion from text, bool, complex or object-valued data."""
    array=np.asarray(value)
    if array.dtype.kind not in 'fiu' or not np.isfinite(array).all():
        raise ValueError('strict finite real numeric operands required')
    return np.asarray(array,float)


def same_numeric_structure(saved,current):
    if type(saved) is not type(current):
        raise ValueError('receipt numeric type or shape changed')
    if isinstance(saved,list):
        if len(saved)!=len(current):
            raise ValueError('receipt numeric type or shape changed')
        for a,b in zip(saved,current):same_numeric_structure(a,b)
    elif type(saved) not in (float,int):
        raise ValueError('receipt numeric type or shape changed')


def exact_record_equal(first,second):
    """Receipt equality includes JSON types, so bool never aliases integer1."""
    if type(first) is not type(second):return False
    if isinstance(first,dict):
        return first.keys()==second.keys() and all(exact_record_equal(first[k],second[k]) for k in first)
    if isinstance(first,list):
        return len(first)==len(second) and all(exact_record_equal(a,b) for a,b in zip(first,second))
    return first==second


def gamma(count):
    eps = np.finfo(float).eps  # Conservative: twice unit roundoff.
    if count*eps >= 1:
        raise ValueError('finite reduction operation count required')
    return count*eps/(1-count*eps)


def weighted_bound(weights, kernels, factor):
    """One product/reduction/outer scale, inflated for its own positive sum."""
    w, k = real_values(weights), real_values(kernels)
    f=real_values(factor)
    if f.shape!=():raise ValueError('scalar real weighted factor required')
    if w.ndim != 1 or k.shape != (len(w), 4) or not np.isfinite(w).all() or not np.isfinite(k).all():
        raise ValueError('finite four-component weighted operands required')
    g = gamma(4*len(w)+8)
    return g*abs(factor)*np.sum(abs(w[:, None]*k), axis=0)*(1+g)


def sum_bound(values):
    rows = [real_values(v) for v in values]
    return gamma(len(rows)+2)*sum((abs(v) for v in rows), np.zeros(4))


def require_values(saved, current, one_replay_bound):
    a, b = real_values(saved), real_values(current)
    limit = 2*real_values(one_replay_bound)
    if np.any(limit<0):raise ValueError('nonnegative arithmetic envelope required')
    limit = np.where(limit, np.nextafter(limit, np.inf), 0.)
    if (a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all()
            or np.any(abs(a-b) > limit)):
        raise ValueError('weighted receipt exceeds replay arithmetic envelope')


def require_record(saved, current, bounds):
    """Only named reduction fields may differ; every other field stays exact."""
    adjusted = deepcopy(current)
    for path, bound in bounds.items():
        old, new, target = saved, current, adjusted
        for key in path[:-1]:
            old, new, target = old[key], new[key], target[key]
        key = path[-1]
        same_numeric_structure(old[key],new[key])
        require_values(old[key], new[key], bound)
        target[key] = old[key]
    if not exact_record_equal(adjusted,saved):
        raise ValueError('weighted receipt differs outside declared reduction fields')
    return saved


def action_bound(value, input_bound):
    from .nsc_incoming_joint_constraints import source_action_gradient
    from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
    jets = incoming_cauchy_jets()
    matrix = np.column_stack([source_action_gradient(jets, row) for row in np.eye(4)])
    return abs(matrix)@input_bound + gamma(32)*(abs(matrix)@abs(np.asarray(value)))


def germ_receipt_bounds(record, payload, joint, baseline):
    """Propagate saved reference insertion sums through the fixed germ algebra."""
    eps=np.finfo(float).eps
    def reference(reference):
        values=[];errors=[]
        for sample in reference['quadrature_samples']:
            w=np.asarray(sample['weights']);insertions=np.asarray(sample['insertions'])
            factor=sample['multiplicity']/(2*np.pi)
            for row in insertions:
                values.append(factor*np.einsum('n,nv->v',w,row))
                errors.append(weighted_bound(w,row,factor))
        return sum(values,np.zeros(4)),sum(errors,np.zeros(4))+sum_bound(values)
    controls=payload['controls'];u=controls['u'];mixed=controls['mixed']
    u16,e16=reference(u['references']['16']);u24,e24=reference(u['references']['24']);mix,emix=reference(mixed['reference'])
    base=np.asarray(baseline['local_action_gradient'])
    def response(local,ref,error):
        return error[:2]+gamma(4)*(abs(np.asarray(local))+abs(base)+abs(ref[:2]))
    eu=response(u['local']['local_action_gradient'],u24,e24)
    ev=response(joint['local_probe']['local_action_gradient'],np.asarray(joint['reference_controls']['n24']['reference_action_gradient_change']),np.zeros(4))
    em=response(mixed['local']['local_action_gradient'],mix,emix)
    changes={name:{(row['field'],row['normal_order'],row['spatial_order']):row['delta'] for row in control['probe']} for name,control in controls.items()}
    pu=changes['u']['r',3,0];mu=changes['mixed']['r',3,0];mv=changes['mixed']['r',1,1]
    coeff=record['numerical_coefficients'];cu,cv2,cv=(coeff[k] for k in ('c_u','c_v2','c_v'))
    ecu=eu[0]/abs(pu)+gamma(3)*abs(cu)
    ecv2=ev[0]/.01**2+gamma(4)*abs(cv2)
    ecv=ev[1]/.01+gamma(3)*abs(cv)
    epred=np.array([abs(mu)*ecu+mv**2*ecv2+gamma(8)*(abs(cu*mu)+abs(cv2*mv**2)),abs(mv)*ecv+gamma(3)*abs(cv*mv)])
    ediff=em+epred+eps*(abs(np.asarray(record['mixed_control_response']))+abs(np.asarray(record['mixed_control_prediction'])))
    bounds={
        ('u_probe_response',):eu,('v_probe_response',):ev,('mixed_control_response',):em,
        ('mixed_control_prediction',):epred,('mixed_control_residual',):ediff,
        ('numerical_coefficients','c_u'):ecu,('numerical_coefficients','c_v2'):ecv2,('numerical_coefficients','c_v'):ecv,
        ('baseline_Jacobian_from_structure',):np.diag([ecu,ecv]),
        ('numerical_Jacobian_determinant',):abs(cv)*ecu+abs(cu)*ecv+ecu*ecv+eps*abs(cu*cv),
        ('residuals','reference_node_replay'):max(float(np.max(e16)),float(np.max(e24)),float(np.max(emix)))
            +eps*max(float(np.max(abs(u16))),float(np.max(abs(u24))),float(np.max(abs(mix)))),
        ('residuals','reference16_24_difference'):float(np.max(e16+e24+eps*(abs(u16)+abs(u24)))),
        ('residuals','homogeneous_shift_zero_identity'):eu[1],
        ('residuals','independent_mixed_control'):float(np.max(ediff)),
    }
    if max(record['residuals'].values()) > record['numerical_tolerance']:
        raise ValueError('original germ numerical control failed')
    return bounds


def receipt_bounds(record, arrays, metadata):
    """Propagate exact saved recipes' weighted sums, additions and subtractions."""
    bounds = {}
    eps = np.finfo(float).eps
    def weighted(prefix, factor, kernel='kernels'):
        return weighted_bound(arrays[prefix+'weights'], arrays[prefix+kernel], factor)
    def add(path, error):
        bounds[tuple(path.split('/'))] = error
    def subtotal(prefixes, factor, kernel='kernels'):
        error = sum((weighted(p, factor, kernel) for p in prefixes), np.zeros(4))
        values = [factor*np.einsum('n,nv->v', arrays[p+'weights'], arrays[p+kernel]) for p in prefixes]
        return error+sum_bound(values)
    def binary(a, b, error_a, error_b):
        return error_a+error_b+eps*(abs(np.asarray(a))+abs(np.asarray(b)))

    schema = record['schema']
    if schema == 'NSC-INCOMING-LOW-HIGH-SOURCE-v1':
        f = record['group_factor_per_sign']; report = metadata['report']
        old = {n: subtotal([f'{s}/archived{n}/' for s in (1,-1)], f) for n in (24,32)}
        new = sum_bound([report['reports'][f'{s}/new64/dps70']['physical_vacuum_source'] for s in (1,-1)])
        for n in (24,32):
            add(f'archived_LOW{n}_source', old[n])
            add(f'new_vacuum_minus_archived_LOW{n}', binary(record['new_order24_vacuum_source'], record[f'archived_LOW{n}_source'], new, old[n]))
        add('archived_LOW32_minus_LOW24', binary(record['archived_LOW32_source'],record['archived_LOW24_source'],old[32],old[24]))
        add('new_order24_vacuum_source', new)
        delta = bounds[('new_vacuum_minus_archived_LOW32',)]
        add('additive_action_gradient_correction', action_bound(record['new_vacuum_minus_archived_LOW32'], delta))
        contraction = max(float(np.max(weighted(prefix+'/', f, 'vacuum_kernels')))
                          for prefix in report['reports'] if '/' in prefix)
        add('residuals/quadrature_contraction_replay', contraction)
    elif schema == 'NSC-INCOMING-MIDDLE-ORDER-CORRECTION-v1':
        f = record['group_factor_per_angular_sign']; report = metadata['report']
        prefixes = [f'{s}/archived16/' for s in (1,-1)]
        old = subtotal(prefixes, f); vacuum = subtotal(prefixes,f,'vacuum_kernels'); thermal = subtotal(prefixes,f,'thermal_kernels')
        original = sum_bound([report['reports'][f'{s}/original24/dps70']['physical_correction'] for s in (1,-1)])
        new = sum_bound([report['reports'][f'{s}/refined48/dps70']['physical_correction'] for s in (1,-1)])
        for field, error in (('old16_middle_source',old),('old16_middle_vacuum',vacuum),('unchanged_archived_thermal_insertion',thermal),('original24_node_vacuum_correction',original),('refined48_node_vacuum_correction',new)):
            add(field,error)
        add('additively_corrected_middle_approximation',binary(record['old16_middle_source'],record['refined48_node_vacuum_correction'],old,new))
        add('source_action_gradient_correction',action_bound(record['refined48_node_vacuum_correction'],new))
        contraction=max(float(np.max(weighted(prefix+'/',f,'correction_kernels'))) for prefix in report['reports'] if '/' in prefix)
        add('residuals/quadrature_contraction_replay',contraction)
    elif schema == 'NSC-INCOMING-GROUP12-ORDER-CORRECTION-v1':
        f = metadata['factor_per_sign']; prefixes=[f'{s}/original16/' for s in (1,-1)]
        old=subtotal(prefixes,f); thermal=subtotal(prefixes,f,'thermal_kernels')
        new=subtotal([f'{s}/refined48/dps70/' for s in (1,-1)],f)
        add('old16_middle_approximant',old);add('retained_thermal_insertion',thermal);add('vacuum_order24_minus16',new)
        add('corrected_middle_approximant',binary(record['old16_middle_approximant'],record['vacuum_order24_minus16'],old,new))
        add('action_gradient_correction',action_bound(record['vacuum_order24_minus16'],new))
        add('maxima/artifact_contraction',max(float(np.max(weighted(prefix+'/',f))) for prefix in metadata['reports'] if '/' in prefix))
    elif schema in ('NSC-INCOMING-WINDOW-REUSE-v1','NSC-INCOMING-GROUP32-SOURCE-v1'):
        for name, window in record['windows'].items():
            f=window['factor_per_sign']; report=metadata['windows'][name]['report']; signs=report['signs'] if 'signs' in report else (1,-1)
            prefixes=[f'{name}/{s}/archived/' for s in signs]
            old=subtotal(prefixes,f)
            new=sum_bound([report['reports'][f'{s}/refined64/dps70']['physical_value'] for s in signs])
            original=sum_bound([report['reports'][f'{s}/original/dps70']['physical_value'] for s in signs])
            if window['kind']=='low_vacuum':
                delta=binary(window['new_numerical_value'],window['old_source'],new,old);updated=new
            else:
                delta=new;updated=binary(window['old_source'],window['new_numerical_value'],old,new)
                add(f'windows/{name}/retained_archived_thermal',subtotal(prefixes,f,'thermal_kernels'))
            for field,error in (('old_source',old),('new_numerical_value',new),('original_grid_numerical_value',original),('explicit_source_delta',delta),('updated_source_approximant',updated)):
                add(f'windows/{name}/{field}',error)
            add(f'windows/{name}/additive_action_gradient',action_bound(window['explicit_source_delta'],delta))
            contraction=max(float(np.max(weighted(f'{name}/{prefix}/',f))) for prefix in report['reports'] if '/' in prefix)
            add(f'windows/{name}/residuals/quadrature_replay',contraction)
        if schema == 'NSC-INCOMING-GROUP32-SOURCE-v1':
            aggregate={}
            for field in ('old_source','explicit_source_delta','updated_source_approximant'):
                aggregate[field]=sum((bounds[('windows',name,field)] for name in record['windows']),np.zeros(4)) +sum_bound([w[field] for w in record['windows'].values()])
                add('aggregate/'+field,aggregate[field])
            add('aggregate/additive_action_gradient',action_bound(record['aggregate']['explicit_source_delta'],aggregate['explicit_source_delta']))
            composition=sum(aggregate.values(),np.zeros(4))+gamma(4)*sum((abs(np.asarray(record['aggregate'][field])) for field in aggregate),np.zeros(4))
            add('aggregate/source_composition_residual',float(np.max(composition)))
    else:
        raise ValueError('explicit weighted receipt recipe required')
    return bounds
