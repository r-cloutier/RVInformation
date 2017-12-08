from RVInformation import *

def assign_observable_flags(self):
    '''
    Assign a flag to each TESS planet highlighting whether or not the planet 
    is reasonably observable with RVs. The conditions for each are defined 
    in the included functions.
    observable_flag == 0: planet is observable in RV
    '''
    # initialize all planets as observable
    self.observable_flag = np.zeros(self.nstars)

    # 
    self.observable_flag = _texp_is_short_enough(self.observable_flag)


def _texp_is_short_enough()
