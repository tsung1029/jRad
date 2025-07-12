def calc_radiation(input_data, track, detector_ene=None, detector_spec=None, myid=0, comm=None):
    import time
    from .track_utils import (
        setup_track,
        read_track,
        track_selection,
        calc_beta,
        calc_beta_dot,
        cleanup_track
    )
    from .calc_spec import calc_spec
    from .calc_energy import calc_ene

    timesum = 0.0
    num_invalid_tracks = 0
    last_track = False

    for p in range(input_data.npart):
        if p == input_data.npart - 1:
            last_track = True

        invalid_track = setup_track(input_data.trackfile, track, input_data.ndimtrack, comm)

        if invalid_track:
            num_invalid_tracks += 1
            if p < input_data.npart - 1:
                continue
            else:
                break

        read_track(input_data.trackfile, track, p, comm)

        track_selection(
            track,
            input_data.track_select_type,
            input_data.nbegin,
            input_data.nend,
            input_data.nrange,
            input_data.x1min,
            input_data.x1max,
            input_data.x1range,
            input_data.tmin,
            input_data.tmax,
            input_data.trange,
            input_data.enemin,
            comm
        )

        if input_data.coherent or not input_data.m_weight:
            track.charge = 1.0

        calc_beta(track)
        calc_beta_dot(track)

        if input_data.diag_type in ("standard", "farfield", "farfieldEndPoints"):
            start_time = time.time()
            calc_spec(track, detector_spec, input_data, last_track)
            end_time = time.time()
            timesum += (end_time - start_time)

        elif input_data.diag_type == "energy":
            calc_ene(track, detector_ene, input_data.diag_type, myid)

        else:
            if myid == 0:
                print("diag_type not available !!")

        cleanup_track(track)

    if myid == 0 and num_invalid_tracks > 0:
        print(f"* warning * Found {num_invalid_tracks} invalid tracks and ignored them\n")