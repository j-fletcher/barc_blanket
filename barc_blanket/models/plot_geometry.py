import matplotlib.pyplot as plt
from barc_model_simple_toroidal import make_model

fig, ax = plt.subplots(1,2, figsize=[12,6])

plot_config = {"section_angle": 180}

model = make_model(plot_config)
universe = model.geometry.root_universe

universe.plot(pixels=500000,
               color_by='cell',
               axes=ax[0]
               )

universe.plot(width=(500.0, 500.0), 
              origin=(680.0, 0.0, 0.1), 
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

plt.savefig('geometry.png')