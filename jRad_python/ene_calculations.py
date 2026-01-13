import numpy as np


class AuxVecsPow:
    def __init__(self, size):
        self.a = np.zeros(size)
        self.b = np.zeros(size)
        self.c = np.zeros(size)
        self.d = np.zeros(size)
        self.n1MinusBeta1xbetaDot2 = np.zeros(size)
        self.n1MinusBeta1xbetaDot3 = np.zeros(size)
        self.n2MinusBeta2xbetaDot1 = np.zeros(size)
        self.n2MinusBeta2xbetaDot3 = np.zeros(size)
        self.n3MinusBeta3xbetaDot1 = np.zeros(size)
        self.n3MinusBeta3xbetaDot2 = np.zeros(size)
        self.n1 = np.zeros(size)
        self.n2 = np.zeros(size)
        self.n3 = np.zeros(size)
        self.r = np.zeros(size)
        self.oneOverR = np.zeros(size)
        self.temp = np.zeros(size)



def calc_ene(track, detector, diag_type: str, myid: int):
    """
    Dispatcher function for calculating energy diagnostics.
    
    Parameters:
    - track: t_track object
    - detector: t_detector_ene object (modified in place)
    - diag_type: str, diagnostic type ('energy', 'power', 'angular_power_dist')
    - myid: int, process ID or identifier
    """

    if diag_type == 'energy':
        calc_radiated_energy(track, detector, myid)

    elif diag_type == 'power':
        print("ERROR: power diagnostic not available")
        raise RuntimeError("Power diagnostic not implemented.")

    elif diag_type == 'angular_power_dist':
        print("ERROR: angular power distribution not available")
        raise RuntimeError("Angular power distribution diagnostic not implemented.")

    else:
        print("ERROR: diag_type not available")
        raise ValueError(f"Diagnostic type '{diag_type}' is not supported.")

def calc_radiated_energy(track, detector, myid):
    ndimtrack = track.ndimtrack
    tracksize = track.tracksize
    ncells1 = detector.ncells1
    ncells2 = detector.ncells2
    x1detmin = detector.x1detmin
    x1detmax = detector.x1detmax
    x2detmin = detector.x2detmin
    x2detmax = detector.x2detmax
    detector_axis = detector.detector_axis.strip()
    x0, y0, z0 = detector.x0, detector.y0, detector.z0
    mynx = detector.xaxis.mynx

    vecs = AuxVecsPow(tracksize)



def calc_radiated_energy(track, detector, myid):
    ndimtrack = track.ndimtrack
    tracksize = track.tracksize
    ncells1 = detector.ncells1
    ncells2 = detector.ncells2
    x1detmin = detector.x1detmin
    x1detmax = detector.x1detmax
    x2detmin = detector.x2detmin
    x2detmax = detector.x2detmax
    detector_axis = detector.detector_axis.strip()
    mynx = detector.xaxis.mynx

    x0 = detector.x0
    y0 = detector.y0
    z0 = detector.z0

    # Allocate arrays
    vecs = {}
    vecs['a'] = np.zeros(tracksize)
    vecs['b'] = np.zeros(tracksize)
    vecs['c'] = np.zeros(tracksize)
    vecs['d'] = np.zeros(tracksize)

    vecs['n1MinusBeta1xbetaDot2'] = np.zeros(tracksize)
    vecs['n1MinusBeta1xbetaDot3'] = np.zeros(tracksize)
    vecs['n2MinusBeta2xbetaDot1'] = np.zeros(tracksize)
    vecs['n2MinusBeta2xbetaDot3'] = np.zeros(tracksize)
    vecs['n3MinusBeta3xbetaDot1'] = np.zeros(tracksize)
    vecs['n3MinusBeta3xbetaDot2'] = np.zeros(tracksize)

    vecs['n1'] = np.zeros(tracksize)
    vecs['n2'] = np.zeros(tracksize)
    vecs['n3'] = np.zeros(tracksize)

    vecs['r'] = np.zeros(tracksize)
    vecs['oneOverR'] = np.zeros(tracksize)
    vecs['temp'] = np.zeros(tracksize)

    if detector.ndim == 1:
        print("Not implemented yet in 1D")

    elif detector.ndim == 2:
        dx1det = (x1detmax - x1detmin) / ncells1
        dx2det = (x2detmax - x2detmin) / ncells2

        for indcell2 in range(ncells2):
            for indcell1 in range(mynx):
                if detector_axis == "x1x2":
                    xi = detector.xaxis.xaxislocal[indcell1]
                    xj = x2detmin + (indcell2 + 0.5) * dx2det
                    xk = z0
                elif detector_axis == "x1x3":
                    xi = detector.xaxis.xaxislocal[indcell1]
                    xj = y0
                    xk = x2detmin + (indcell2 + 0.5) * dx2det
                elif detector_axis == "x2x3":
                    xi = x0
                    xj = detector.xaxis.xaxislocal[indcell1]
                    xk = x2detmin + (indcell2 + 0.5) * dx2det
                else:
                    raise ValueError("Invalid detector axis")

                if ndimtrack == 1:
                    vecs['r'] = np.sqrt((track.x1 - xi)**2 + xj**2 + xk**2)
                    vecs['n1'] = (xi - track.x1) / vecs['r']
                    vecs['n2'] = xj / vecs['r']
                    vecs['n3'] = xk / vecs['r']
                elif ndimtrack == 2:
                    vecs['r'] = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + xk**2)
                    vecs['n1'] = (xi - track.x1) / vecs['r']
                    vecs['n2'] = (xj - track.x2) / vecs['r']
                    vecs['n3'] = xk / vecs['r']
                elif ndimtrack == 3:
                    vecs['r'] = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + (track.x3 - xk)**2)
                    vecs['n1'] = (xi - track.x1) / vecs['r']
                    vecs['n2'] = (xj - track.x2) / vecs['r']
                    vecs['n3'] = (xk - track.x3) / vecs['r']
                else:
                    raise ValueError("ndimtrack must be 1, 2, or 3")

                calc_radiated_energy_cell(track, detector, vecs, indcell1, indcell2)

    else:
        print("Error: not implemented")

def calc_radiated_energy_cell(track, detector, vecs, indcell1, indcell2):
    import numpy as np

    dx1det = (detector.x1detmax - detector.x1detmin) / detector.ncells1
    dx2det = (detector.x2detmax - detector.x2detmin) / detector.ncells2

    tracksize = track.tracksize

    vecs.n1MinusBeta1xbetaDot2 = (vecs.n1 - track.beta1) * track.betaDot2
    vecs.n1MinusBeta1xbetaDot3 = (vecs.n1 - track.beta1) * track.betaDot3
    vecs.n2MinusBeta2xbetaDot1 = (vecs.n2 - track.beta2) * track.betaDot1
    vecs.n2MinusBeta2xbetaDot3 = (vecs.n2 - track.beta2) * track.betaDot3
    vecs.n3MinusBeta3xbetaDot1 = (vecs.n3 - track.beta3) * track.betaDot1
    vecs.n3MinusBeta3xbetaDot2 = (vecs.n3 - track.beta3) * track.betaDot2

    # Calculate power per solid angle dP/dOmega
    vecs.a =  vecs.n2 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
              vecs.n3 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1)
    vecs.b = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
              vecs.n3 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)
    vecs.c = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1) - \
              vecs.n2 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)

    vecs.d = 1.0 - (vecs.n1 * track.beta1 + vecs.n2 * track.beta2 + vecs.n3 * track.beta3)

    power_density = (vecs.a**2 + vecs.b**2 + vecs.c**2) * track.dt[0] * (dx1det * dx2det) / ((vecs.r**2) * (vecs.d**5))
    detector.pow2d[indcell1, indcell2] += abs(track.charge) * np.sum(power_density)

