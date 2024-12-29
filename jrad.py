# jrad_def.py
# ================
# Define the basic data structure for the jRad library


import numpy as np
import re
import copy as cp
from fractions import Fraction as frac


class Input_class:
    def __init__(self, wmin, wmax, wpoints,waxistype,min_dec,num_dec,ppdec, trackfile
        , ndimtrack, track_select_type, nbegin, nend, nrange, x1min, x1max, x1range,
        tmin, tmax, trange, enemin, npart, m_weight, ndim, diag_type, endpoints, detector_axis
        , ncells1, ncells2, x0, y0, z0, x1detmin, x1detmax, x2detmin, x2detmax, emissivity
        , coherent):
        self.wmin = wmin
        self.wmax = wmax
        self.wpoints = wpoints
        self.waxistype = waxistype
        self.min_dec = min_dec
        self.num_dec = num_dec
        self.ppdec = ppdec
        self.trackfile = trackfile
        self.ndimtrack = ndimtrack
        self.track_select_type = track_select_type
        self.nbegin = nbegin
        self.nend = nend
        self.nrange = nrange
        self.x1min = x1min
        self.x1max = x1max
        self.x1range = x1range
        self.tmin = tmin
        self.tmax = tmax
        self.trange = trange
        self.enemin = enemin
        self.npart = npart
        self.m_weight = m_weight
        self.ndim = ndim
        self.diag_type = diag_type
        self.endpoints = endpoints
        self.detector_axis = detector_axis
        self.ncells1 = ncells1
        self.ncells2 = ncells2
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.x1detmin = x1detmin
        self.x1detmax = x1detmax
        self.x2detmin = x2detmin
        self.x2detmax = x2detmax
        self.emissivity = emissivity
        self.coherent = coherent


def read_input_file(filename=None):

    import f90nml

    if filename is None:
        filename='input'

    # waxis_parameters
    wmin = -1.0
    wmax = -1.0
    wpoints=-1
    waxistype = '-'
    min_dec=0
    num_dec=0
    ppdec=-1

    # track_parameters
    trackfile = '-'
    ndimtrack=-1
    track_select_type = 'none'
    nbegin=-1
    nend = -1
    # nrange(1:2) = 0
    nrange = [0,0]
    x1min = 0.0
    x1max = 0.0
    x1range = [0.0, 0.0]
    tmin = 0.0
    tmax = 0.0
    trange = [ 0.0, 0.0],
    enemin = 0.0
    npart=-1
    m_weight = False

    # detector_parameters
    ndim = -1
    diag_type='-'
    endpoints=True
    detector_axis = ''
    ncells1 = -1
    ncells2=-1
    x0 = -1.0
    y0 = -1.0
    z0 = -1.0
    x1detmin = 0.0
    x1detmax = 0.0
    x2detmin = 0.0
    x2detmax = 0.0
    emissivity = 'd2W/dwdO'
    coherent = False

    input_file_data=f90nml.read(filename)

    if (input_file_data['waxis_parameters']['wmin'] is not None):
        wmin = input_file_data['waxis_parameters']['wmin']
    if (input_file_data['waxis_parameters']['wmax'] is not None):
        wmax = input_file_data['waxis_parameters']['wmax']
    if (input_file_data['waxis_parameters']['wpoints'] is not None):
        wpoints = input_file_data['waxis_parameters']['wpoints']
    if (input_file_data['waxis_parameters']['waxistype'] is not None):
        waxistype = input_file_data['waxis_parameters']['waxistype']
    # if (input_file_data['waxis_parameters']['min_dec'] is not None):
    if('min_dec' in input_file_data['waxis_parameters']):
        min_dec = input_file_data['waxis_parameters']['min_dec']
    # if (input_file_data['waxis_parameters']['num_dec'] is not None):
    if('num_dec' in input_file_data['waxis_parameters']):
        num_dec = input_file_data['waxis_parameters']['num_dec']
    # if (input_file_data['waxis_parameters']['ppdec'] is not None):
    if('ppdec' in input_file_data['waxis_parameters']):
        ppdec = input_file_data['waxis_parameters']['ppdec']

    if (input_file_data['track_parameters']['trackfile'] is not None):
        trackfile = input_file_data['track_parameters']['trackfile']
    if (input_file_data['track_parameters']['ndimtrack'] is not None):
        ndimtrack = input_file_data['track_parameters']['ndimtrack']
    if (input_file_data['track_parameters']['track_select_type'] is not None):
        track_select_type = input_file_data['track_parameters']['track_select_type']
    if (input_file_data['track_parameters']['nbegin'] is not None):
        nbegin = input_file_data['track_parameters']['nbegin']
    if (input_file_data['track_parameters']['nend'] is not None):
        nend = input_file_data['track_parameters']['nend']
    if (input_file_data['track_parameters']['nrange'] is not None):
        nrange = input_file_data['track_parameters']['nrange']
    if (input_file_data['track_parameters']['x1min'] is not None):
        x1min = input_file_data['track_parameters']['x1min']
    if (input_file_data['track_parameters']['x1max'] is not None):
        x1max = input_file_data['track_parameters']['x1max']
    if (input_file_data['track_parameters']['x1range'] is not None):
        x1range = input_file_data['track_parameters']['x1range']
    if (input_file_data['track_parameters']['tmin'] is not None):
        tmin = input_file_data['track_parameters']['tmin']
    if (input_file_data['track_parameters']['tmax'] is not None):
        tmax = input_file_data['track_parameters']['tmax']
    if (input_file_data['track_parameters']['trange'] is not None):
        trange = input_file_data['track_parameters']['trange']
    if (input_file_data['track_parameters']['enemin'] is not None):
        enemin = input_file_data['track_parameters']['enemin']
    if (input_file_data['track_parameters']['npart'] is not None):
        npart = input_file_data['track_parameters']['npart']
    if (input_file_data['track_parameters']['m_weight'] is not None):
        m_weight = input_file_data['track_parameters']['m_weight']


    if (input_file_data['detector_parameters']['ndim'] is not None):
        ndim = input_file_data['detector_parameters']['ndim']
    if (input_file_data['detector_parameters']['diag_type'] is not None):
        diag_type = input_file_data['detector_parameters']['diag_type']
    if (input_file_data['detector_parameters']['endpoints'] is not None):
        endpoints = input_file_data['detector_parameters']['endpoints']
    if (input_file_data['detector_parameters']['detector_axis'] is not None):
        detector_axis = input_file_data['detector_parameters']['detector_axis']
    if (input_file_data['detector_parameters']['ncells1'] is not None):
        ncells1 = input_file_data['detector_parameters']['ncells1']
    if (input_file_data['detector_parameters']['ncells2'] is not None):
        ncells2 = input_file_data['detector_parameters']['ncells2']
    # if (input_file_data['detector_parameters']['x0'] is not None):
    if('x0' in input_file_data['detector_parameters']):
        x0 = input_file_data['detector_parameters']['x0']
    # if (input_file_data['detector_parameters']['y0'] is not None):
    if('y0' in input_file_data['detector_parameters']):
        y0 = input_file_data['detector_parameters']['y0']
    # if (input_file_data['detector_parameters']['z0'] is not None):
    if('z0' in input_file_data['detector_parameters']):
        z0 = input_file_data['detector_parameters']['z0']
    # if (input_file_data['detector_parameters']['x1detmin'] is not None):
    if('x1detmin' in input_file_data['detector_parameters']):
        x1detmin = input_file_data['detector_parameters']['x1detmin']
    # if (input_file_data['detector_parameters']['x1detmax'] is not None):
    if('x1detmax' in input_file_data['detector_parameters']):
        x1detmax = input_file_data['detector_parameters']['x1detmax']
    # if (input_file_data['detector_parameters']['x2detmin'] is not None):
    if('x2detmin' in input_file_data['detector_parameters']):
        x2detmin = input_file_data['detector_parameters']['x2detmin']
    # if (input_file_data['detector_parameters']['x2detmax'] is not None):
    if('x2detmax' in input_file_data['detector_parameters']):
        x2detmax = input_file_data['detector_parameters']['x2detmax']
    # if (input_file_data['detector_parameters']['emissivity'] is not None):
    if('emissivity' in input_file_data['detector_parameters']):
        emissivity = input_file_data['detector_parameters']['emissivity']
    # if (input_file_data['detector_parameters']['coherent'] is not None):
    if('coherent' in input_file_data['detector_parameters']):
        coherent = input_file_data['detector_parameters']['coherent']

# here we initialize the variable from the input file
    input_var = Input_class(
        wmin = wmin,
        wmax = wmax,
        wpoints = wpoints,
        waxistype = waxistype,
        min_dec = min_dec,
        num_dec = num_dec,
        ppdec = ppdec,
        trackfile = trackfile,
        ndimtrack = ndimtrack,
        track_select_type = track_select_type,
        nbegin = nbegin,
        nend = nend,
        nrange = nrange,
        x1min = x1min,
        x1max = x1max,
        x1range = x1range,
        tmin = tmin,
        tmax = tmax,
        trange = trange,
        enemin = enemin,
        npart = npart,
        m_weight = m_weight,
        ndim = ndim,
        diag_type = diag_type,
        endpoints = endpoints,
        detector_axis = detector_axis,
        ncells1 = ncells1,
        ncells2 = ncells2,
        x0 = x0,
        y0 = y0,
        z0 = z0,
        x1detmin = x1detmin,
        x1detmax = x1detmax,
        x2detmin = x2detmin,
        x2detmax = x2detmax,
        emissivity = emissivity,
        coherent = coherent)

    return input_var




class Track_class:
    def __init__(self, tracksize=tracksize,ndimtrack=ndimtrack):
        self.tracksize=tracksize
        self.ndimtrack=ndimtrack
        self.t=np.zeros(tracksize,dtype = np.float64,'F')
        self.x1=np.zeros(tracksize,dtype = np.float64,'F')
        self.x2=np.zeros(tracksize,dtype = np.float64,'F')
        self.x3=np.zeros(tracksize,dtype = np.float64,'F')
        self.p1=np.zeros(tracksize,dtype = np.float64,'F')
        self.p2=np.zeros(tracksize,dtype = np.float64,'F')
        self.p3=np.zeros(tracksize,dtype = np.float64,'F')
        self.t=np.zeros(tracksize,dtype = np.float64,'F')
        self.t=np.zeros(tracksize,dtype = np.float64,'F')
        self.t=np.zeros(tracksize,dtype = np.float64,'F')
        self.beta1=np.zeros(tracksize,dtype=np.float64,'F')
        self.beta2=np.zeros(tracksize,dtype=np.float64,'F')
        self.beta3=np.zeros(tracksize,dtype=np.float64,'F')
        self.betaDot1=np.zeros(tracksize,dtype=np.float64,'F')
        self.betaDot2=np.zeros(tracksize,dtype=np.float64,'F')
        self.betaDot3=np.zeros(tracksize,dtype=np.float64,'F')

        self.g=np.zeros(tracksize,dtype=np.float64,'F')

class waxis_class:
    def __init__(self,input_var,detector_var):
        # we are trying to replicate the functions of MAKE_WAXIS in jrad.f90
        # 


class t_detector_ene_class:

    def__init__(self,detector_var,npart,m_weight):

    




