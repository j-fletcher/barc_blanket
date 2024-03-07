import openmc


def evaluate_metric(model:openmc.Model, metric):
    """ Evaluate the metric for the given model

    Parameters:
    ----------
    model : openmc.Model
        The model to evaluate the metric for
    metric : str
        The name of the metric to evaluate

    Returns:
    -------
    metric_val: float
        The value of the metric calculated for the model
    """

    if metric == "tbr":
        metric_val = tritium_breeding_ratio(model)
    #elif metric == "some_other_arbitrary_metric":
    #    metric_val = whatever_function(model)
    else:
        raise ValueError(f"Invalid metric: {metric}")

    return metric_val

# THIS IS JUST A PROOF OF CONCEPT
# THIS FUNCTION DOES NOT PRODUCE CORRECT OUTPUT
def tritium_breeding_ratio(model:openmc.Model):
    """ THIS IS JUST A PROOF OF CONCEPT
    THIS FUNCTION DOES NOT PRODUCE CORRECT OUTPUT
    
    Calculate the tritium breeding ratio for the given model

    Parameters:
    ----------
    model : openmc.Model
        The model to evaluate the metric for

    Returns:
    -------
    tbr: float
        The tritium breeding ratio for the model
    """

    # Run the model
    model.run()
    final_statepoint = openmc.StatePoint("statepoint.50.h5")

    # Get tally results
    # TODO: do this programmatically instead of hardcoded here
    tallies = final_statepoint.tallies
    tbr_tally_id = 2
    tally_result = tallies[tbr_tally_id].mean[0][0][0]
    
    # TODO do some volume weighting or whatever to get an actual TBR
    return tally_result