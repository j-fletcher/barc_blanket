import matplotlib.pyplot as plt

def plot_geometry(model):
    fig, ax = plt.subplots(1,2, figsize=[12,6])

    universe = model.geometry.root_universe

    universe.plot(pixels=500000,
                  width=(800,400),
                  origin=(480, 150, 0),
                  color_by='cell',
                axes=ax[0]
                )

    universe.plot(width=(800.0, 800.0), 
                origin=(480.0, 0.0, 0.1), 
                basis='xz',
                pixels=500000,
                color_by='cell',
                axes=ax[1]
                )

    ax[0].set_title('Top View')
    ax[0].set_xlabel('x(cm)')
    ax[0].set_ylabel('y(cm)')
    ax[1].set_title('Azimuthal View')
    ax[1].set_xlabel('x(cm)')
    ax[1].set_ylabel('z(cm)')
    ax[1].legend(loc='upper right')
    fig.tight_layout()

    return fig