"""Certified reduction of an already enclosed axial Fourier profile table.

Discarded coefficient balls are added to the continuous derivative tail.
This changes the numerical representation, never the physical profile.
"""
from flint import arb,ctx


def compress_profiles(profiles, period_length, retained_index, *, bits=120):
    if isinstance(retained_index,bool) or not isinstance(retained_index,int) or retained_index<1:
        raise ValueError('positive retained Fourier index required')
    with ctx.workprec(bits):
        length=arb(period_length)
        if not length>0:raise ValueError('positive period required')
        omega=2*arb.pi()/length
        result={}
        for key,profile in profiles.items():
            coefficients=profile['coefficients']
            old=(len(coefficients)-1)//2
            if len(coefficients)!=2*old+1 or retained_index>old:
                raise ValueError('compression must stay inside the existing complete band')
            tail=[arb(v) for v in profile['tail']]
            if any(not v.is_finite() or not v>=0 for v in tail):
                raise ValueError('finite nonnegative prior tails required')
            for j,value in enumerate(coefficients):
                index=j-old
                if abs(index)>retained_index:
                    magnitude=value.abs_upper()
                    rate=omega*abs(index)
                    for order in range(len(tail)):
                        tail[order] += rate**order*magnitude
            result[key]={**profile,
                'coefficients':list(coefficients[old-retained_index:old+retained_index+1]),
                'tail':tuple(v.upper() for v in tail),
                'retained_index':retained_index,'original_retained_index':old,
                'discarded_coefficients_bounded':True}
        return result
