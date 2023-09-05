import numpy as np
from tqdm import tqdm


def load_event_field(event, kind, paths, grid, MAS="PCS", in_rsp=True,
                     smooth_scales=[0.0, 1.0, 2.0, 3.0, 4.0]):
    """
    Load a field evaluated at the locations of a given event.

    Parameters
    ----------
    event : str
        Event name.
    kind : str
        Field type.
    paths : gwlss.Paths
        Paths object.
    grid : int
        Grid size.
    MAS : str, optional
         Mass-assignment scheme.
    in_rsp : bool, optional
        Whether the field is in redshift space.
    smooth_scales : list, optional
        List of smoothing scales that must be present in the pre-computed
        field.

    Returns
    -------
    3-dimensional array of shape `(nsims, nsamples, nsmooth)`
    """
    nsims = [7444 + n * 24 for n in range(101)]

    for n, nsim in enumerate(tqdm(nsims)):

        fpath = paths.evaluated_field(event, kind, nsim, grid, MAS, in_rsp)

        with np.load(fpath) as f:
            field_val = f["val"]
            field_smooth_scales = f["smooth_scales"]

        if n == 0:
            for smooth_scale in smooth_scales:
                if smooth_scale not in field_smooth_scales:
                    raise ValueError(f"Smooth scale {smooth_scale} not found.")

            ks = np.where(np.isin(field_smooth_scales, smooth_scales))[0]

            shape = (len(nsims), len(field_val), len(smooth_scales))
            out = np.zeros(shape, dtype=field_val.dtype)

        out[n, ...] = field_val[:, ks]

    return out
