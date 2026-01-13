#
# spec calculations.py, 
# chatgpt translated from spec_calculations.f90
#

@dataclass
class TAuxVecs:
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    d: Optional[np.ndarray] = None

    n1: np.ndarray = field(default_factory=lambda: np.array([]))
    n2: np.ndarray = field(default_factory=lambda: np.array([]))
    n3: np.ndarray = field(default_factory=lambda: np.array([]))
    ndotbeta: np.ndarray = field(default_factory=lambda: np.array([]))
    ndotr: Optional[np.ndarray] = None
    r: np.ndarray = field(default_factory=lambda: np.array([]))
    faccorr: Optional[np.ndarray] = None

    n1MinusBeta1xbetaDot2: Optional[np.ndarray] = None
    n1MinusBeta1xbetaDot3: Optional[np.ndarray] = None
    n2MinusBeta2xbetaDot1: Optional[np.ndarray] = None
    n2MinusBeta2xbetaDot3: Optional[np.ndarray] = None
    n3MinusBeta3xbetaDot1: Optional[np.ndarray] = None
    n3MinusBeta3xbetaDot2: Optional[np.ndarray] = None

    temp1cRe: np.ndarray = field(default_factory=lambda: np.array([]))
    temp1cIm: np.ndarray = field(default_factory=lambda: np.array([]))
    temp2cRe: np.ndarray = field(default_factory=lambda: np.array([]))
    temp2cIm: np.ndarray = field(default_factory=lambda: np.array([]))
    temp3cRe: np.ndarray = field(default_factory=lambda: np.array([]))
    temp3cIm: np.ndarray = field(default_factory=lambda: np.array([]))

    exppartRe: np.ndarray = field(default_factory=lambda: np.array([]))
    exppartIm: np.ndarray = field(default_factory=lambda: np.array([]))
    commonpart1Re: np.ndarray = field(default_factory=lambda: np.array([]))
    commonpart2Re: np.ndarray = field(default_factory=lambda: np.array([]))
    commonpart2Im: np.ndarray = field(default_factory=lambda: np.array([]))

    temp: np.ndarray = field(default_factory=lambda: np.array([]))


def calc_spec(track, detector, input_obj, last_track):
    diag_type = input_obj.diag_type.lower()

    if diag_type == 'standard':
        calc_spec_standard(track, detector, input_obj, last_track)
    elif diag_type in ('farfield', 'farfieldendpoints'):
        calc_spec_farfield(track, detector, input_obj, last_track)
    else:
        raise ValueError("ERROR: diag_type not available")

def calc_spec_standard(track, detector, input_obj, last_track):
    ndimtrack = track.ndimtrack
    tracksize = track.tracksize
    ncells1 = detector.ncells1
    ncells2 = detector.ncells2
    x1detmin = detector.x1detmin
    x1detmax = detector.x1detmax
    x2detmin = detector.x2detmin
    x2detmax = detector.x2detmax
    detector_axis = detector.detector_axis.strip()

    # allocate auxiliary arrays
    vecs = t_auxvecs()
    vecs.allocate(tracksize)

    # zero out required arrays
    vecs.zero(['exppartRe', 'exppartIm',
               'temp1cRe', 'temp1cIm',
               'temp2cRe', 'temp2cIm',
               'temp3cRe', 'temp3cIm',
               'commonpart1Re', 'commonpart2Re', 'commonpart2Im'])

    x0 = detector.x0
    y0 = detector.y0
    z0 = detector.z0

    if detector.ndim == 1:
        xi, xj, xk = x0, y0, z0
        if ndimtrack == 1:
            vecs.r = np.sqrt((track.x1 - xi)**2 + xj**2 + xk**2)
            vecs.n1 = (xi - track.x1) / vecs.r
            vecs.n2 = xj / vecs.r
            vecs.n3 = xk / vecs.r
        elif ndimtrack == 2:
            vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + xk**2)
            vecs.n1 = (xi - track.x1) / vecs.r
            vecs.n2 = (xj - track.x2) / vecs.r
            vecs.n3 = xk / vecs.r
        elif ndimtrack == 3:
            vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + (track.x3 - xk)**2)
            vecs.n1 = (xi - track.x1) / vecs.r
            vecs.n2 = (xj - track.x2) / vecs.r
            vecs.n3 = (xk - track.x3) / vecs.r
        else:
            raise ValueError("Error: ndimtrack must be 1, 2, or 3")

        calc_spec_standard_cell(track, detector, input_obj, vecs, 0, 0, last_track)

    elif detector.ndim == 2:
        dx1det = (x1detmax - x1detmin) / ncells1
        for indcell1 in range(1, ncells1 + 1):
            if detector_axis == "x1":
                xi = x1detmin + (indcell1 - 0.5) * dx1det
                xj, xk = y0, z0
            elif detector_axis == "x2":
                xi = x0
                xj = x1detmin + (indcell1 - 0.5) * dx1det
                xk = z0
            elif detector_axis == "x3":
                xi, xj = x0, y0
                xk = x1detmin + (indcell1 - 0.5) * dx1det

            if ndimtrack == 1:
                vecs.r = np.sqrt((track.x1 - xi)**2 + xj**2 + xk**2)
                vecs.n1 = (xi - track.x1) / vecs.r
                vecs.n2 = xj / vecs.r
                vecs.n3 = xk / vecs.r
            elif ndimtrack == 2:
                vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + xk**2)
                vecs.n1 = (xi - track.x1) / vecs.r
                vecs.n2 = (xj - track.x2) / vecs.r
                vecs.n3 = xk / vecs.r
            elif ndimtrack == 3:
                vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + (track.x3 - xk)**2)
                vecs.n1 = (xi - track.x1) / vecs.r
                vecs.n2 = (xj - track.x2) / vecs.r
                vecs.n3 = (xk - track.x3) / vecs.r
            else:
                raise ValueError("Error: ndimtrack must be 1, 2, or 3")

            calc_spec_standard_cell(track, detector, input_obj, vecs, indcell1, 0, last_track)

    elif detector.ndim == 3:
        dx1det = (x1detmax - x1detmin) / ncells1
        dx2det = (x2detmax - x2detmin) / ncells2

        for indcell2 in range(1, ncells2 + 1):
            for indcell1 in range(1, ncells1 + 1):
                if detector_axis == "x1x2":
                    xi = x1detmin + (indcell1 - 0.5) * dx1det
                    xj = x2detmin + (indcell2 - 0.5) * dx2det
                    xk = z0
                elif detector_axis == "x1x3":
                    xi = x1detmin + (indcell1 - 0.5) * dx1det
                    xj = y0
                    xk = x2detmin + (indcell2 - 0.5) * dx2det
                elif detector_axis == "x2x3":
                    xi = x0
                    xj = x1detmin + (indcell1 - 0.5) * dx1det
                    xk = x2detmin + (indcell2 - 0.5) * dx2det

                if ndimtrack == 1:
                    vecs.r = np.sqrt((track.x1 - xi)**2 + xj**2 + xk**2)
                    vecs.n1 = (xi - track.x1) / vecs.r
                    vecs.n2 = xj / vecs.r
                    vecs.n3 = xk / vecs.r
                elif ndimtrack == 2:
                    vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + xk**2)
                    vecs.n1 = (xi - track.x1) / vecs.r
                    vecs.n2 = (xj - track.x2) / vecs.r
                    vecs.n3 = xk / vecs.r
                elif ndimtrack == 3:
                    vecs.r = np.sqrt((track.x1 - xi)**2 + (track.x2 - xj)**2 + (track.x3 - xk)**2)
                    vecs.n1 = (xi - track.x1) / vecs.r
                    vecs.n2 = (xj - track.x2) / vecs.r
                    vecs.n3 = (xk - track.x3) / vecs.r
                else:
                    raise ValueError("Error: ndimtrack must be 1, 2, or 3")

                calc_spec_standard_cell(track, detector, input_obj, vecs, indcell1, indcell2, last_track)

    # deallocate (implicitly handled in Python, or define vecs.deallocate() if needed)
    vecs.deallocate()


def calc_spec_standard_cell(track, detector, input_obj, vecs, indcell1, indcell2, last_track):
    tracksize = track.tracksize

    # Vector precomputations
    vecs.n1MinusBeta1xbetaDot2 = aMinusBResTimesC(vecs.n1, track.beta1, track.betaDot2, tracksize)
    vecs.n1MinusBeta1xbetaDot3 = aMinusBResTimesC(vecs.n1, track.beta1, track.betaDot3, tracksize)
    vecs.n2MinusBeta2xbetaDot1 = aMinusBResTimesC(vecs.n2, track.beta2, track.betaDot1, tracksize)
    vecs.n2MinusBeta2xbetaDot3 = aMinusBResTimesC(vecs.n2, track.beta2, track.betaDot3, tracksize)
    vecs.n3MinusBeta3xbetaDot1 = aMinusBResTimesC(vecs.n3, track.beta3, track.betaDot1, tracksize)
    vecs.n3MinusBeta3xbetaDot2 = aMinusBResTimesC(vecs.n3, track.beta3, track.betaDot2, tracksize)

    vecs.ndotbeta = vecs.n1 * track.beta1 + vecs.n2 * track.beta2 + vecs.n3 * track.beta3

    vecs.a = vecs.n2 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
             vecs.n3 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1)

    vecs.b = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
              vecs.n3 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)

    vecs.c = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1) - \
              vecs.n2 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)

    vecs.d = 1.0 - vecs.ndotbeta

    # First part of common product
    vecs.commonpart1Re = calcCommonpart1RePart1(vecs.d, track.dt, tracksize)

    # Begin frequency loop
    wwBegin = 2 if detector.waxis.waxislocal[0] == 0.0 else 1

    for ww in range(wwBegin - 1, detector.waxis.mynw):  # Python uses 0-based indexing
        w = detector.waxis.waxislocal[ww]
        vecs.temp = w * (track.t + vecs.r)

        vecs.exppartRe = np.cos(vecs.temp)
        vecs.exppartIm = np.sin(vecs.temp)

        vecs.commonpart2Re = vecs.commonpart1Re * vecs.exppartRe
        vecs.commonpart2Im = vecs.commonpart1Re * vecs.exppartIm

        # Component-wise multiplication and trapezoidal correction
        for comp, name in zip([vecs.a, vecs.b, vecs.c], ['temp1c', 'temp2c', 'temp3c']):
            re = vecs.commonpart2Re * comp
            im = vecs.commonpart2Im * comp

            re[0] *= 0.5
            re[-1] *= 0.5
            im[0] *= 0.5
            im[-1] *= 0.5

            setattr(vecs, f'{name}Re', re)
            setattr(vecs, f'{name}Im', im)

        # Apply emissivity option
        if input_obj.emissivity == "d2W/dwdS":
            for name in ['temp1c', 'temp2c', 'temp3c']:
                setattr(vecs, f'{name}Re', getattr(vecs, f'{name}Re') / vecs.r)
                setattr(vecs, f'{name}Im', getattr(vecs, f'{name}Im') / vecs.r)

        # Integrate via sum
        temp1csumRe = np.sum(vecs.temp1cRe)
        temp1csumIm = np.sum(vecs.temp1cIm)
        temp2csumRe = np.sum(vecs.temp2cRe)
        temp2csumIm = np.sum(vecs.temp2cIm)
        temp3csumRe = np.sum(vecs.temp3cRe)
        temp3csumIm = np.sum(vecs.temp3cIm)

        # Fill spectrum
        ndim = detector.ndim
        coherent = input_obj.coherent

        if coherent:
            if ndim == 1:
                detector.spec_1d_c1_re[ww] += temp1csumRe
                detector.spec_1d_c1_im[ww] += temp1csumIm
                detector.spec_1d_c2_re[ww] += temp2csumRe
                detector.spec_1d_c2_im[ww] += temp2csumIm
                detector.spec_1d_c3_re[ww] += temp3csumRe
                detector.spec_1d_c3_im[ww] += temp3csumIm

                if last_track:
                    detector.spec1d[ww] += (
                        detector.spec_1d_c1_re[ww] ** 2 + detector.spec_1d_c1_im[ww] ** 2 +
                        detector.spec_1d_c2_re[ww] ** 2 + detector.spec_1d_c2_im[ww] ** 2 +
                        detector.spec_1d_c3_re[ww] ** 2 + detector.spec_1d_c3_im[ww] ** 2
                    ) / 4.0

            elif ndim == 2:
                reim = [
                    (detector.spec_2d_c1_re, temp1csumRe), (detector.spec_2d_c1_im, temp1csumIm),
                    (detector.spec_2d_c2_re, temp2csumRe), (detector.spec_2d_c2_im, temp2csumIm),
                    (detector.spec_2d_c3_re, temp3csumRe), (detector.spec_2d_c3_im, temp3csumIm),
                ]
                for arr, val in reim:
                    arr[ww, indcell1 - 1] += val

                if last_track:
                    c1 = detector.spec_2d_c1_re[ww, indcell1 - 1] ** 2 + detector.spec_2d_c1_im[ww, indcell1 - 1] ** 2
                    c2 = detector.spec_2d_c2_re[ww, indcell1 - 1] ** 2 + detector.spec_2d_c2_im[ww, indcell1 - 1] ** 2
                    c3 = detector.spec_2d_c3_re[ww, indcell1 - 1] ** 2 + detector.spec_2d_c3_im[ww, indcell1 - 1] ** 2
                    detector.spec2d[ww, indcell1 - 1] += (c1 + c2 + c3) / 4.0

            elif ndim == 3:
                i, j = indcell1 - 1, indcell2 - 1
                reim = [
                    (detector.spec_3d_c1_re, temp1csumRe), (detector.spec_3d_c1_im, temp1csumIm),
                    (detector.spec_3d_c2_re, temp2csumRe), (detector.spec_3d_c2_im, temp2csumIm),
                    (detector.spec_3d_c3_re, temp3csumRe), (detector.spec_3d_c3_im, temp3csumIm),
                ]
                for arr, val in reim:
                    arr[ww, i, j] += val

                if last_track:
                    c1 = detector.spec_3d_c1_re[ww, i, j] ** 2 + detector.spec_3d_c1_im[ww, i, j] ** 2
                    c2 = detector.spec_3d_c2_re[ww, i, j] ** 2 + detector.spec_3d_c2_im[ww, i, j] ** 2
                    c3 = detector.spec_3d_c3_re[ww, i, j] ** 2 + detector.spec_3d_c3_im[ww, i, j] ** 2
                    detector.spec3d[ww, i, j] += (c1 + c2 + c3) / 4.0

        else:
            # Incoherent case
            spectrum_val = abs(track.charge) * (
                temp1csumRe**2 + temp1csumIm**2 +
                temp2csumRe**2 + temp2csumIm**2 +
                temp3csumRe**2 + temp3csumIm**2
            ) / 4.0

            if ndim == 1:
                detector.spec1d[ww] += spectrum_val
            elif ndim == 2:
                detector.spec2d[ww, indcell1 - 1] += spectrum_val
            elif ndim == 3:
                detector.spec3d[ww, indcell1 - 1, indcell2 - 1] += spectrum_val

# 🧱 You Will Need to Implement:
# aMinusBResTimesC(a, b, c, size) → element-wise a - b * c

# calcCommonpart1RePart1(d, dt, size) → function returning array of size size

# vecs: should support fields as NumPy arrays and allow field assignment as shown

def calc_spec_standard_cell_ver2(track, detector, input, vecs, indcell1, indcell2, last_track):
    tracksize = track.tracksize
    waxislocal = detector.waxis.waxislocal
    mynw = detector.waxis.mynw

    # Vector calculations (skipping aMinusBResTimesC for now)
    vecs.ndotbeta = vecs.n1 * track.beta1 + vecs.n2 * track.beta2 + vecs.n3 * track.beta3
    vecs.a = vecs.n2 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
             vecs.n3 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1)
    vecs.b = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot2 - vecs.n2MinusBeta2xbetaDot1) + \
              vecs.n3 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)
    vecs.c = -vecs.n1 * (vecs.n1MinusBeta1xbetaDot3 - vecs.n3MinusBeta3xbetaDot1) - \
              vecs.n2 * (vecs.n2MinusBeta2xbetaDot3 - vecs.n3MinusBeta3xbetaDot2)
    vecs.d = 1.0 - vecs.ndotbeta

    vecs.commonpart1Re = vecs.d * track.dt  # calcCommonpart1RePart1 approximation

    wwBegin = 1 if waxislocal[0] != 0.0 else 2
    for ww in range(wwBegin - 1, mynw):
        w = waxislocal[ww]
        vecs.temp = w * (track.t + vecs.r)
        vecs.exppartRe = np.cos(vecs.temp)
        vecs.exppartIm = np.sin(vecs.temp)

        vecs.commonpart2Re = vecs.commonpart1Re * vecs.exppartRe
        vecs.commonpart2Im = vecs.commonpart1Re * vecs.exppartIm

        vecs.temp1cRe = vecs.commonpart2Re * vecs.a
        vecs.temp1cIm = vecs.commonpart2Im * vecs.a
        vecs.temp2cRe = vecs.commonpart2Re * vecs.b
        vecs.temp2cIm = vecs.commonpart2Im * vecs.b
        vecs.temp3cRe = vecs.commonpart2Re * vecs.c
        vecs.temp3cIm = vecs.commonpart2Im * vecs.c

        for temp in [vecs.temp1cRe, vecs.temp1cIm, vecs.temp2cRe, vecs.temp2cIm, vecs.temp3cRe, vecs.temp3cIm]:
            temp[0] *= 0.5
            temp[-1] *= 0.5

        if input.emissivity == "d2W/dwdS":
            vecs.temp1cRe /= vecs.r
            vecs.temp1cIm /= vecs.r
            vecs.temp2cRe /= vecs.r
            vecs.temp2cIm /= vecs.r
            vecs.temp3cRe /= vecs.r
            vecs.temp3cIm /= vecs.r

        temp1re = trapezoidal_sum(vecs.temp1cRe)
        temp1im = trapezoidal_sum(vecs.temp1cIm)
        temp2re = trapezoidal_sum(vecs.temp2cRe)
        temp2im = trapezoidal_sum(vecs.temp2cIm)
        temp3re = trapezoidal_sum(vecs.temp3cRe)
        temp3im = trapezoidal_sum(vecs.temp3cIm)

        if input.coherent:
            if detector.ndim == 1:
                accumulate_coherent_1d(detector, ww, temp1re, temp1im, temp2re, temp2im, temp3re, temp3im, last_track)
            # 2D and 3D cases to be implemented similarly
        else:
            if detector.ndim == 1:
                accumulate_incoherent_1d(detector, ww, temp1re, temp1im, temp2re, temp2im, temp3re, temp3im, track.charge)
            # 2D and 3D cases to be implemented similarly



def zero_vector(vec: np.ndarray) -> None:
    vec.fill(0.0)


def a_minus_b_res_times_c(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return a * (1.0 - b * c)

def calc_commonpart1_re_part1(d: np.ndarray, dt: float) -> np.ndarray:
    return dt / d

def trapezoidal_sum(y: np.ndarray) -> float:
    y_copy = y.copy()
    y_copy[0] *= 0.5
    y_copy[-1] *= 0.5
    return np.sum(y_copy)

def accumulate_coherent_1d(detector, ww, temp1re, temp1im, temp2re, temp2im, temp3re, temp3im, last_track):
    detector.spec_1d_c1_re[ww] += temp1re
    detector.spec_1d_c1_im[ww] += temp1im
    detector.spec_1d_c2_re[ww] += temp2re
    detector.spec_1d_c2_im[ww] += temp2im
    detector.spec_1d_c3_re[ww] += temp3re
    detector.spec_1d_c3_im[ww] += temp3im

    if last_track:
        detector.spec1d[ww] += (
            (detector.spec_1d_c1_re[ww] ** 2 + detector.spec_1d_c1_im[ww] ** 2) +
            (detector.spec_1d_c2_re[ww] ** 2 + detector.spec_1d_c2_im[ww] ** 2) +
            (detector.spec_1d_c3_re[ww] ** 2 + detector.spec_1d_c3_im[ww] ** 2)
        ) / 4.0

def accumulate_incoherent_1d(detector, ww, temp1re, temp1im, temp2re, temp2im, temp3re, temp3im, charge):
    detector.spec1d[ww] += abs(charge) * (
        (temp1re ** 2 + temp1im ** 2) +
        (temp2re ** 2 + temp2im ** 2) +
        (temp3re ** 2 + temp3im ** 2)
    ) / 4.0

def calc_spec_farfield_cell(track, detector, input_params, vecs, indcell1, indcell2, last_track):
    tracksize = track.tracksize

    # Compute ndotbeta
    vecs.ndotbeta = vecs.n1 * track.beta1 + vecs.n2 * track.beta2 + vecs.n3 * track.beta3

    # Compute vecs.a, b, c
    vecs.a = vecs.n2 * track.betaDot3 - vecs.n3 * track.betaDot2
    vecs.b = vecs.n3 * track.betaDot1 - vecs.n1 * track.betaDot3
    vecs.c = vecs.n1 * track.betaDot2 - vecs.n2 * track.betaDot1

    vecs.d = 1.0 - vecs.ndotbeta

    # commonpart1Re: dt/d * d (vector)
    vecs.commonpart1Re = vecs.faccorr * track.dt / vecs.d

    # Determine start index for frequency loop
    wwBegin = 1 if detector.waxis.waxislocal[0] != 0.0 else 1

    for ww in range(wwBegin, detector.waxis.mynw):
        w = detector.waxis.waxislocal[ww]

        vecs.temp = w * (track.t - vecs.ndotr)

        vecs.exppartRe = np.cos(vecs.temp)
        vecs.exppartIm = np.sin(vecs.temp)

        vecs.commonpart2Re = vecs.commonpart1Re * vecs.exppartRe
        vecs.commonpart2Im = vecs.commonpart1Re * vecs.exppartIm

        vecs.temp1cRe = vecs.commonpart2Re * vecs.a
        vecs.temp1cRe[0] *= 0.5
        vecs.temp1cRe[-1] *= 0.5

        vecs.temp1cIm = vecs.commonpart2Im * vecs.a
        vecs.temp1cIm[0] *= 0.5
        vecs.temp1cIm[-1] *= 0.5

        vecs.temp2cRe = vecs.commonpart2Re * vecs.b
        vecs.temp2cRe[0] *= 0.5
        vecs.temp2cRe[-1] *= 0.5

        vecs.temp2cIm = vecs.commonpart2Im * vecs.b
        vecs.temp2cIm[0] *= 0.5
        vecs.temp2cIm[-1] *= 0.5

        vecs.temp3cRe = vecs.commonpart2Re * vecs.c
        vecs.temp3cRe[0] *= 0.5
        vecs.temp3cRe[-1] *= 0.5

        vecs.temp3cIm = vecs.commonpart2Im * vecs.c
        vecs.temp3cIm[0] *= 0.5
        vecs.temp3cIm[-1] *= 0.5

        temp1csumRe = np.sum(vecs.temp1cRe)
        temp1csumIm = np.sum(vecs.temp1cIm)
        temp2csumRe = np.sum(vecs.temp2cRe)
        temp2csumIm = np.sum(vecs.temp2cIm)
        temp3csumRe = np.sum(vecs.temp3cRe)
        temp3csumIm = np.sum(vecs.temp3cIm)

        if input_params.coherent:
            if detector.ndim == 1:
                detector.spec_1d_c1_re[ww] += temp1csumRe
                detector.spec_1d_c1_im[ww] += temp1csumIm
                detector.spec_1d_c2_re[ww] += temp2csumRe
                detector.spec_1d_c2_im[ww] += temp2csumIm
                detector.spec_1d_c3_re[ww] += temp3csumRe
                detector.spec_1d_c3_im[ww] += temp3csumIm

                if last_track:
                    detector.spec1d[ww] += (
                        detector.spec_1d_c1_re[ww] ** 2 + detector.spec_1d_c1_im[ww] ** 2 +
                        detector.spec_1d_c2_re[ww] ** 2 + detector.spec_1d_c2_im[ww] ** 2 +
                        detector.spec_1d_c3_re[ww] ** 2 + detector.spec_1d_c3_im[ww] ** 2
                    ) / 4.0

            elif detector.ndim == 2:
                detector.spec_2d_c1_re[ww, indcell1] += temp1csumRe
                detector.spec_2d_c1_im[ww, indcell1] += temp1csumIm
                detector.spec_2d_c2_re[ww, indcell1] += temp2csumRe
                detector.spec_2d_c2_im[ww, indcell1] += temp2csumIm
                detector.spec_2d_c3_re[ww, indcell1] += temp3csumRe
                detector.spec_2d_c3_im[ww, indcell1] += temp3csumIm

                if last_track:
                    detector.spec2d[ww, indcell1] += (
                        detector.spec_2d_c1_re[ww, indcell1] ** 2 + detector.spec_2d_c1_im[ww, indcell1] ** 2 +
                        detector.spec_2d_c2_re[ww, indcell1] ** 2 + detector.spec_2d_c2_im[ww, indcell1] ** 2 +
                        detector.spec_2d_c3_re[ww, indcell1] ** 2 + detector.spec_2d_c3_im[ww, indcell1] ** 2
                    ) / 4.0

            elif detector.ndim == 3:
                detector.spec_3d_c1_re[ww, indcell1, indcell2] += temp1csumRe
                detector.spec_3d_c1_im[ww, indcell1, indcell2] += temp1csumIm
                detector.spec_3d_c2_re[ww, indcell1, indcell2] += temp2csumRe
                detector.spec_3d_c2_im[ww, indcell1, indcell2] += temp2csumIm
                detector.spec_3d_c3_re[ww, indcell1, indcell2] += temp3csumRe
                detector.spec_3d_c3_im[ww, indcell1, indcell2] += temp3csumIm

                if last_track:
                    detector.spec3d[ww, indcell1, indcell2] += (
                        detector.spec_3d_c1_re[ww, indcell1, indcell2] ** 2 +
                        detector.spec_3d_c1_im[ww, indcell1, indcell2] ** 2 +
                        detector.spec_3d_c2_re[ww, indcell1, indcell2] ** 2 +
                        detector.spec_3d_c2_im[ww, indcell1, indcell2] ** 2 +
                        detector.spec_3d_c3_re[ww, indcell1, indcell2] ** 2 +
                        detector.spec_3d_c3_im[ww, indcell1, indcell2] ** 2
                    ) / 4.0
        else:
            charge_abs = abs(track.charge)
            if detector.ndim == 1:
                detector.spec1d[ww] += charge_abs * (
                    temp1csumRe ** 2 + temp1csumIm ** 2 +
                    temp2csumRe ** 2 + temp2csumIm ** 2 +
                    temp3csumRe ** 2 + temp3csumIm ** 2
                ) / 4.0
            elif detector.ndim == 2:
                detector.spec2d[ww, indcell1] += charge_abs * (
                    temp1csumRe ** 2 + temp1csumIm ** 2 +
                    temp2csumRe ** 2 + temp2csumIm ** 2 +
                    temp3csumRe ** 2 + temp3csumIm ** 2
                ) / 4.0
            elif detector.ndim == 3:
                detector.spec3d[ww, indcell1, indcell2] += charge_abs * (
                    temp1csumRe ** 2 + temp1csumIm ** 2 +
                    temp2csumRe ** 2 + temp2csumIm ** 2 +
                    temp3csumRe ** 2 + temp3csumIm ** 2
                ) / 4.0

def calc_spec_farfield(track, detector, input_data, last_track):
    from numpy import sqrt, zeros, full

    indcell1 = None
    indcell2 = None
    ndimtrack = track.ndimtrack
    tracksize = track.tracksize
    ncells1 = detector.ncells1
    ncells2 = detector.ncells2
    x1detmin = detector.x1detmin
    x1detmax = detector.x1detmax
    x2detmin = detector.x2detmin
    x2detmax = detector.x2detmax
    detector_axis = detector.detector_axis.strip()
    
    # Allocate and zero vectors
    vecs = AuxVecs(tracksize)
    vecs.faccorr = facQC / track.g

    x0, y0, z0 = detector.x0, detector.y0, detector.z0

    ndim = detector.ndim

    if ndim == 1:
        xi, xj, xk = x0, y0, z0
        vecs.r = sqrt(xi**2 + xj**2 + xk**2)
        vecs.n1 = full(tracksize, xi / vecs.r)
        vecs.n2 = full(tracksize, xj / vecs.r)
        vecs.n3 = full(tracksize, xk / vecs.r)

        if ndimtrack == 1:
            vecs.ndotr = vecs.n1 * track.x1
        elif ndimtrack == 2:
            vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2
        elif ndimtrack == 3:
            vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2 + vecs.n3 * track.x3
        else:
            raise ValueError("ndimtrack must be 1, 2, or 3")

        calc_spec_farfield_cell(track, detector, input_data, vecs, None, None, last_track)

    elif ndim == 2:
        dx1det = (x1detmax - x1detmin) / ncells1

        for indcell1 in range(1, ncells1 + 1):
            if detector_axis == "x1":
                xi = x1detmin + (indcell1 - 0.5) * dx1det
                xj = y0
                xk = z0
            elif detector_axis == "x2":
                xi = x0
                xj = x1detmin + (indcell1 - 0.5) * dx1det
                xk = z0
            elif detector_axis == "x3":
                xi = x0
                xj = y0
                xk = x1detmin + (indcell1 - 0.5) * dx1det

            vecs.r = sqrt(xi**2 + xj**2 + xk**2)
            vecs.n1 = full(tracksize, xi / vecs.r)
            vecs.n2 = full(tracksize, xj / vecs.r)
            vecs.n3 = full(tracksize, xk / vecs.r)

            if ndimtrack == 1:
                vecs.ndotr = vecs.n1 * track.x1
            elif ndimtrack == 2:
                vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2
            elif ndimtrack == 3:
                vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2 + vecs.n3 * track.x3
            else:
                raise ValueError("ndimtrack must be 1, 2, or 3")

            calc_spec_farfield_cell(track, detector, input_data, vecs, indcell1, None, last_track)

    elif ndim == 3:
        dx1det = (x1detmax - x1detmin) / ncells1
        dx2det = (x2detmax - x2detmin) / ncells2

        for indcell2 in range(1, ncells2 + 1):
            for indcell1 in range(1, ncells1 + 1):
                if detector_axis == "x1x2":
                    xi = x1detmin + (indcell1 - 0.5) * dx1det
                    xj = x2detmin + (indcell2 - 0.5) * dx2det
                    xk = z0
                elif detector_axis == "x1x3":
                    xi = x1detmin + (indcell1 - 0.5) * dx1det
                    xj = y0
                    xk = x2detmin + (indcell2 - 0.5) * dx2det
                elif detector_axis == "x2x3":
                    xi = x0
                    xj = x1detmin + (indcell1 - 0.5) * dx1det
                    xk = x2detmin + (indcell2 - 0.5) * dx2det

                vecs.r = sqrt(xi**2 + xj**2 + xk**2)
                vecs.n1 = full(tracksize, xi / vecs.r)
                vecs.n2 = full(tracksize, xj / vecs.r)
                vecs.n3 = full(tracksize, xk / vecs.r)

                if ndimtrack == 1:
                    vecs.ndotr = vecs.n1 * track.x1
                elif ndimtrack == 2:
                    vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2
                elif ndimtrack == 3:
                    vecs.ndotr = vecs.n1 * track.x1 + vecs.n2 * track.x2 + vecs.n3 * track.x3
                else:
                    raise ValueError("ndimtrack must be 1, 2, or 3")

                calc_spec_farfield_cell(track, detector, input_data, vecs, indcell1, indcell2, last_track)

    else:
        raise ValueError("Unsupported detector ndim: should be 1, 2, or 3")
    
