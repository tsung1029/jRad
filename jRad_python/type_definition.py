from dataclasses import dataclass, field
import numpy as np

@dataclass
class TWAxis:
    waxistype: str = ""
    firstwvec: np.ndarray = None  # Assume integer pointer
    waxis: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    waxislocal: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    mynw: int = 0
    nw: int = 0
    dw: float = 0.0
    wmin: float = 0.0
    wmax: float = 0.0
    dimsfwaxis: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=int))
    dims_chunkwaxis: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=int))

@dataclass
class TXAxis:
    firstvec: np.ndarray = None  # Assume integer pointer
    xaxis: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    xaxislocal: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    mynx: int = 0
    nx: int = 0
    dx: float = 0.0
    xmin: float = 0.0
    xmax: float = 0.0
    dimsf_xaxis: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=int))
    dims_chunk_xaxis: np.ndarray = field(default_factory=lambda: np.zeros(1, dtype=int))

@dataclass
class TTrack:
    n: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    t: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    x1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    x2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    x3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    p1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    p2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    p3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    dt: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    ene: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    chargetemp: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    charge: float = 0.0
    tracksize: int = 0
    ndimtrack: int = 0
    dtinitial: float = 0.0
    beta1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    beta2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    beta3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    betaDot1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    betaDot2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    betaDot3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    g: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))

@dataclass
class TDetectorSpec:
    spec1d: np.ndarray = None
    spec2d: np.ndarray = None
    spec3d: np.ndarray = None

    spec_1d_c1_re: np.ndarray = None
    spec_1d_c1_im: np.ndarray = None
    spec_1d_c2_re: np.ndarray = None
    spec_1d_c2_im: np.ndarray = None
    spec_1d_c3_re: np.ndarray = None
    spec_1d_c3_im: np.ndarray = None

    spec_2d_c1_re: np.ndarray = None
    spec_2d_c1_im: np.ndarray = None
    spec_2d_c2_re: np.ndarray = None
    spec_2d_c2_im: np.ndarray = None
    spec_2d_c3_re: np.ndarray = None
    spec_2d_c3_im: np.ndarray = None

    spec_3d_c1_re: np.ndarray = None
    spec_3d_c1_im: np.ndarray = None
    spec_3d_c2_re: np.ndarray = None
    spec_3d_c2_im: np.ndarray = None
    spec_3d_c3_re: np.ndarray = None
    spec_3d_c3_im: np.ndarray = None

    e1_1d_re: np.ndarray = None
    e1_1d_im: np.ndarray = None
    e2_1d_re: np.ndarray = None
    e2_1d_im: np.ndarray = None
    e3_1d_re: np.ndarray = None
    e3_1d_im: np.ndarray = None

    e1_2d_re: np.ndarray = None
    e1_2d_im: np.ndarray = None
    e2_2d_re: np.ndarray = None
    e2_2d_im: np.ndarray = None
    e3_2d_re: np.ndarray = None
    e3_2d_im: np.ndarray = None

    e1_3d_re: np.ndarray = None
    e1_3d_im: np.ndarray = None
    e2_3d_re: np.ndarray = None
    e2_3d_im: np.ndarray = None
    e3_3d_re: np.ndarray = None
    e3_3d_im: np.ndarray = None

    coherent: bool = False
    dimsf: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    dims_chunk: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    ndim: int = 0
    ncells1: int = 0
    ncells2: int = 0
    x0: float = 0.0
    y0: float = 0.0
    z0: float = 0.0
    x1detmin: float = 0.0
    x1detmax: float = 0.0
    x2detmin: float = 0.0
    x2detmax: float = 0.0
    detector_axis: str = ""
    waxis: TWAxis = field(default_factory=TWAxis)

@dataclass
class TDetectorEne:
    pow1d: np.ndarray = None
    pow2d: np.ndarray = None
    dimsf: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    dims_chunk: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    ndim: int = 0
    ncells1: int = 0
    ncells2: int = 0
    x0: float = 0.0
    y0: float = 0.0
    z0: float = 0.0
    x1detmin: float = 0.0
    x1detmax: float = 0.0
    x2detmin: float = 0.0
    x2detmax: float = 0.0
    detector_axis: str = ""
    xaxis: TXAxis = field(default_factory=TXAxis)

@dataclass
class TAuxVecs:
    a: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    b: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    c: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    d: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    ndotr: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1MinusBeta1xbetaDot2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1MinusBeta1xbetaDot3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2MinusBeta2xbetaDot1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2MinusBeta2xbetaDot3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3MinusBeta3xbetaDot1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3MinusBeta3xbetaDot2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    ndotbeta: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    eta: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    r: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    faccorr: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp1cRe: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp1cIm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp2cRe: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp2cIm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp3cRe: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp3cIm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    exppartRe: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    exppartIm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    commonpart1Re: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    commonpart2Re: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    commonpart2Im: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))

@dataclass
class TAuxVecsPow:
    a: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    b: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    c: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    d: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    r: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    oneOverR: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1MinusBeta1xbetaDot2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n1MinusBeta1xbetaDot3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2MinusBeta2xbetaDot1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n2MinusBeta2xbetaDot3: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3MinusBeta3xbetaDot1: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    n3MinusBeta3xbetaDot2: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))
    temp: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float64))

@dataclass
class TInput:
    wmin: float = 0.0
    wmax: float = 0.0
    min_dec: int = 0
    num_dec: int = 0
    ppdec: int = 0
    wpoints: int = 0
    waxistype: str = ""
    trackfile: str = ""
    track_select_type: str = ""
    nbegin: int = 0
    nend: int = 0
    npart: int = 0
    ndimtrack: int = 0
    nrange: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=int))
    x1range: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=np.float64))
    trange: np.ndarray = field(default_factory=lambda: np.zeros(2, dtype=np.float64))
    x1min: float = 0.0
    x1max: float = 0.0
    tmin: float = 0.0
    tmax: float = 0.0
    enemin: float = 0.0
    m_weight: bool = False
    ndim: int = 0
    diag_type: str = ""
    endpoints: bool = False
    detector_axis: str = ""
    ncells1: int = 0
    ncells2: int = 0
    x0: float = 0.0
    y0: float = 0.0
    z0: float = 0.0
    x1detmin: float = 0.0
    x1detmax: float = 0.0
    x2detmin: float = 0.0
    x2detmax: float = 0.0
    emissivity: str = ""
    coherent: bool = False
    filename: str = ""
    parallelIO: bool = False
