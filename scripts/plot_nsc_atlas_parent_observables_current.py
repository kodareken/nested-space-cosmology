"""Current presentation of the pinned observable plate, with parent-count wording.

The historical producer and its numerical/source hashes remain unchanged.
This scoped save hook changes only the identified title before rendering.
"""
from matplotlib.figure import Figure
import plot_nsc_atlas_parent_observables as retained


def main():
    original=Figure.savefig
    def current_view(figure,*args,**kwargs):
        for axis in figure.axes:
            if axis.get_title()=='A  |  Expanding parent, retained 4D exterior':
                axis.set_title('A  |  Expanding spatial parent')
        return original(figure,*args,**kwargs)
    Figure.savefig=current_view
    try:retained.main()
    finally:Figure.savefig=original


if __name__=='__main__':main()
