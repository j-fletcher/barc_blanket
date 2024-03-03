import openmc






def make_model(model_config):

    # Make tallies

    tbr_cell_filter = openmc.CellFilter([fusion_blanket_cell, burner_blanket_cell])

    model=openmc.Model(geometry=geometry,
                       settings=settings,
                       tallie=tallies)

    return model