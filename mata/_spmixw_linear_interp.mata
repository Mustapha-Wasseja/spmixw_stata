*! version 0.1.0  2026-05-08  spmixw: 1D linear interpolation
version 13.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------
// _spmixw_linear_interp()
//
// Piecewise-linear interpolation of (x_in, y_in) onto x_out, with edge
// clamping (values outside [min(x_in), max(x_in)] return the nearest
// endpoint). Mirrors R's `approx(..., rule = 2)`.
//
// x_in must be sorted ascending.
// ----------------------------------------------------------------------
real colvector _spmixw_linear_interp(real colvector x_in,
                                     real colvector y_in,
                                     real colvector x_out)
{
    real scalar    nx, no, i, j
    real scalar    t
    real colvector y_out

    nx = rows(x_in)
    no = rows(x_out)
    y_out = J(no, 1, 0)

    // Single-pass walk: x_out is assumed sorted ascending and the coarse
    // grid x_in is sorted. We carry j across iterations so each interp
    // is amortised O(1) and the whole call is O(nx + no), not O(nx*no).
    j = 1
    for (i = 1; i <= no; i++) {
        if (x_out[i] <= x_in[1]) {
            y_out[i] = y_in[1]
        }
        else if (x_out[i] >= x_in[nx]) {
            y_out[i] = y_in[nx]
        }
        else {
            while (j < nx - 1 & x_in[j+1] < x_out[i]) {
                j = j + 1
            }
            t        = (x_out[i] - x_in[j]) / (x_in[j+1] - x_in[j])
            y_out[i] = y_in[j] + t * (y_in[j+1] - y_in[j])
        }
    }

    return(y_out)
}

end
