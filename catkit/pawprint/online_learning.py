import numpy as np
import matplotlib.pyplot as plt
import sklearn.gaussian_process as gp


class OnlineLeaener():
    """ online_learner

    A class for incremental gaussian process learning. This will
    learn from the available labelled data and from the computed 
    acquisition function it will output the point that needs to 
    be evaluated next. Based on the chosen acquisition function, 
    it will foucs on exploration and exploitation of the machine 
    learning space. 

    """

    def __init__(self,
            kernel=None,
            n_batches=10,
            batch_policy='equal_spacing',
            acquisition_function='prob_improv', 
            )
