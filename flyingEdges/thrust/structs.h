#ifndef STRUCTS_STUFF
#define STRUCTS_STUFF

#include "config.h"

struct abc_transform
  : public thrust::unary_function<
        int,
        tuple<
            tuple<
                scalar_t, scalar_t, scalar_t,
                scalar_t, scalar_t, scalar_t>,
            tuple<
                scalar_t, scalar_t, scalar_t,
                scalar_t, scalar_t, scalar_t>,
            tuple<
                scalar_t, scalar_t, scalar_t,
                scalar_t, scalar_t, scalar_t>,
            uchar,
            int> >
{
    abc_transform(
        int const& nx, int const& ny, int const& nz,
        int const& spacing_x, int const& spacing_y, int const& spacing_z,
        int const& zeropos_x, int const& zeropos_y, int const& zeropos_z,
        scalar_t const& isoval,
        const_pointer<scalar_t> data)
      : nx(nx), ny(ny), nz(nz),
        spacing_x(spacing_x), spacing_y(spacing_y), spacing_z(spacing_z),
        zeropos_x(zeropos_x), zeropos_y(zeropos_y), zeropos_z(zeropos_z),
        isoval(isoval),
        data(data)
    {}

    int const nx, ny, nz;
    int const spacing_x, spacing_y, spacing_z;
    int const zeropos_x, zeropos_y, zeropos_z;
    scalar_t const isoval;
    const_pointer<scalar_t> data;

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    tuple<
        tuple<
            scalar_t, scalar_t, scalar_t,    // points  ax, ay, az
            scalar_t, scalar_t, scalar_t>,   // normals ax, ay, az
        tuple<
            scalar_t, scalar_t, scalar_t,    // points  bx, by, bz
            scalar_t, scalar_t, scalar_t>,   // normals bx, by, bz
        tuple<
            scalar_t, scalar_t, scalar_t,    // points  cx, cy, cz
            scalar_t, scalar_t, scalar_t>,   // normals cx, cy, cz
        uchar,                               // cube_id
        int>                                 // nt
    operator()(
        int const& idx) const
    {
        int const k = idx / (nx*ny);
        int const j = (idx / nx) % ny;
        int const i = idx % nx;

        scalar_t vals[8];
        vals[0] = get_safe(i,     j,     k);
        vals[1] = get_safe(i + 1, j,     k);
        vals[2] = get_safe(i + 1, j + 1, k);
        vals[3] = get_safe(i,     j + 1, k);
        vals[4] = get_safe(i,     j,     k + 1);
        vals[5] = get_safe(i + 1, j,     k + 1);
        vals[6] = get_safe(i + 1, j + 1, k + 1);
        vals[7] = get_safe(i,     j + 1, k + 1);

        uchar cube_id = calc_cube_id(vals);
        int num_tri = mctable::numTris[cube_id];

        if(i == nx-1 || j == ny-1 || k == nz-1)
        {
            // don't set cube_id to 0, might still need to calc a, b or c
            num_tri = 0;
        }

        bool const* is_cut = mctable::isCut[cube_id];

        bool calc_a = (i != nx-1) && is_cut[0];
        bool calc_b = (j != ny-1) && is_cut[3];
        bool calc_c = (k != nz-1) && is_cut[8];

        scalar_t points[3*3];  // x,y,z coords for a,b,c
        scalar_t normals[3*3]; // x,y,z coords for a,b,c

        calc_points(
            points,
            vals,
            i, j, k,
            calc_a, calc_b, calc_c);
        calc_normals(
            normals,
            vals,
            i, j, k,
            calc_a, calc_b, calc_c);

        // get min pos for possible points and set that minus 1 as the
        // not set value.
        scalar_t default_value = zeropos_x;
        if(zeropos_y < default_value)
            default_value = zeropos_y;
        if(zeropos_z < default_value)
            default_value = zeropos_z;

        default_value = default_value - 1;
        auto default_tuple =
            make_tuple(
                default_value, default_value, default_value,
                default_value, default_value, default_value);

        return make_tuple(
            calc_a
              ? make_tuple(
                    points[0],  points[1],  points[2],
                    normals[0], normals[1], normals[2])
              : default_tuple,
            calc_b
              ? make_tuple(
                    points[3],  points[4],  points[5],
                    normals[3], normals[4], normals[5])
              : default_tuple,
            calc_c
              ? make_tuple(
                    points[6],  points[7],  points[8],
                    normals[6], normals[7], normals[8])
              : default_tuple,
            cube_id,
            num_tri);
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    uchar
    calc_cube_id(scalar_t* vals) const
    {
        uchar case_id = 0;

        if(vals[0] >= isoval) case_id |= 1;
        if(vals[1] >= isoval) case_id |= 2;
        if(vals[2] >= isoval) case_id |= 4;
        if(vals[3] >= isoval) case_id |= 8;
        if(vals[4] >= isoval) case_id |= 16;
        if(vals[5] >= isoval) case_id |= 32;
        if(vals[6] >= isoval) case_id |= 64;
        if(vals[7] >= isoval) case_id |= 128;

        return case_id;
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    void
    calc_points(
        scalar_t* pts,
        scalar_t* vals,
        int const& i, int const& j, int const& k,
        bool const& calc_a, bool const& calc_b, bool const& calc_c) const
    {
        scalar_t xpos = zeropos_x + i * spacing_x;
        scalar_t ypos = zeropos_y + j * spacing_y;
        scalar_t zpos = zeropos_z + k * spacing_z;

        if(calc_a)
        {
            scalar_t weight = (isoval - vals[0]) / (vals[1] - vals[0]);
            pts[0] = interpolate(weight, xpos, xpos + spacing_x);
            pts[1] = ypos;
            pts[2] = zpos;
        }

        if(calc_b)
        {
            scalar_t weight = (isoval - vals[0]) / (vals[3] - vals[0]);
            pts[3] = xpos;
            pts[4] = interpolate(weight, ypos, ypos + spacing_y);
            pts[5] = zpos;
        }

        if(calc_c)
        {
            scalar_t weight = (isoval - vals[0]) / (vals[4] - vals[0]);
            pts[6] = xpos;
            pts[7] = ypos;
            pts[8] = interpolate(weight, zpos, zpos + spacing_z);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    void
    calc_normals(
        scalar_t* nrs,
        scalar_t* vals,
        int const& i, int const& j, int const& k,
        bool const& calc_a, bool const& calc_b, bool const& calc_c) const
    {
        scalar_t g0x, g0y, g0z;
        compute_gradient(
            i, j, k,
            g0x, g0y, g0z);

        if(calc_a)
        {
            scalar_t g1x, g1y, g1z;
            compute_gradient(
                i+1, j, k,
                g1x, g1y, g1z);

            scalar_t weight = (isoval - vals[0]) / (vals[1] - vals[0]);
            nrs[0] = interpolate(weight, g0x, g1x);
            nrs[1] = interpolate(weight, g0y, g1y);
            nrs[2] = interpolate(weight, g0z, g1z);
        }

        if(calc_b)
        {
            scalar_t g3x, g3y, g3z;
            compute_gradient(
                i, j+1, k,
                g3x, g3y, g3z);

            scalar_t weight = (isoval - vals[0]) / (vals[3] - vals[0]);
            nrs[3] = interpolate(weight, g0x, g3x);
            nrs[4] = interpolate(weight, g0y, g3y);
            nrs[5] = interpolate(weight, g0z, g3z);
        }

        if(calc_c)
        {
            scalar_t g4x, g4y, g4z;
            compute_gradient(
                i, j, k+1,
                g4x, g4y, g4z);

            scalar_t weight = (isoval - vals[0]) / (vals[4] - vals[0]);
            nrs[6] = interpolate(weight, g0x, g4x);
            nrs[7] = interpolate(weight, g0y, g4y);
            nrs[8] = interpolate(weight, g0z, g4z);
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    void
    compute_gradient(
        int const& i, int const& j, int const& k,
        scalar_t& gx, scalar_t& gy, scalar_t& gz) const
    {
        scalar_t x[3][2];
        scalar_t run[3];

        if (i == 0)
        {
            x[0][0] = get(i+1, j, k);
            x[0][1] = get(i, j, k);
            run[0] = spacing_x;
        }
        else if (i == (nx - 1))
        {
            x[0][0] = get(i, j, k);
            x[0][1] = get(i-1, j, k);
            run[0] = spacing_x;
        }
        else
        {
            x[0][0] = get(i+1, j, k);
            x[0][1] = get(i-1, j, k);
            run[0] = 2 * spacing_x;
        }

        if (j == 0)
        {
            x[1][0] = get(i, j+1, k);
            x[1][1] = get(i, j, k);
            run[1] = spacing_y;
        }
        else if (j == (ny - 1))
        {
            x[1][0] = get(i, j, k);
            x[1][1] = get(i, j-1, k);
            run[1] = spacing_y;
        }
        else
        {
            x[1][0] = get(i, j+1, k);
            x[1][1] = get(i, j-1, k);
            run[1] = 2 * spacing_y;
        }

        if (k == 0)
        {
            x[2][0] = get(i, j, k+1);
            x[2][1] = get(i, j, k);
            run[2] = spacing_z;
        }
        else if (k == (nz - 1))
        {
            x[2][0] = get(i, j, k);
            x[2][1] = get(i, j, k-1);
            run[2] = spacing_z;
        }
        else
        {
            x[2][0] = get(i, j, k+1);
            x[2][1] = get(i, j, k-1);
            run[2] = 2 * spacing_z;
        }

        gx = (x[0][1] - x[0][0]) / run[0];
        gy = (x[1][1] - x[1][0]) / run[1];
        gz = (x[2][1] - x[2][0]) / run[2];
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    scalar_t
    interpolate(
        scalar_t const& weight,
        scalar_t const& lhs,
        scalar_t const& rhs) const
    {
        return lhs + (weight * (rhs - lhs));
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    scalar_t
    get_safe(int const& i, int const& j, int const& k) const
    {
        if(i == nx || j == ny || k == nz)
            return *(data + idxer(i-1, j-1, k-1)); // arbitrary value, doesn't matter
        return *(data + idxer(i, j, k));
    }

    __platform__
    scalar_t
    get(int const& i, int const& j, int const& k) const
    {
        return *(data + idxer(i, j, k));
    }

    ///////////////////////////////////////////////////////////////////////////
    __platform__
    int
    idxer(int const& i, int const& j, int const& k) const
    {
        return k*nx*ny + j*nx + i;
    }

};

struct neq
  : public unary_function<
        scalar_t,
        int>
{
    neq(scalar_t const& default_value)
      : default_value(default_value)
    {}

    scalar_t default_value;

    __platform__
    int
    operator()(scalar_t const& val) const
    {
        return val != default_value;
    }
};

struct set_points_and_normals
  : public unary_function<
        tuple<
            tuple<
                int, int, int>,
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  ax, ay, az
                scalar_t, scalar_t, scalar_t>,   // normals ax, ay, az
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  bx, by, bz
                scalar_t, scalar_t, scalar_t>,   // normals bx, by, bz
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  cx, cy, cz
                scalar_t, scalar_t, scalar_t> >, // normals cx, cy, cz
            void>
{
    set_points_and_normals(
        pointer<scalar_t> pts_x,
        pointer<scalar_t> pts_y,
        pointer<scalar_t> pts_z,
        pointer<scalar_t> nrs_x,
        pointer<scalar_t> nrs_y,
        pointer<scalar_t> nrs_z,
        scalar_t const& default_value,
        int const& offset)
      : pts_x(pts_x), pts_y(pts_y), pts_z(pts_z),
        nrs_x(nrs_x), nrs_y(nrs_y), nrs_z(nrs_z),
        default_value(default_value),
        offset(offset)
    {}

    pointer<scalar_t> pts_x, pts_y, pts_z;
    pointer<scalar_t> nrs_x, nrs_y, nrs_z;
    scalar_t const default_value;
    int const offset;

    __device__ __host__
    void
    operator()(
        tuple<
            tuple<
                int, int, int>, // a0, b0, c0
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  ax, ay, az
                scalar_t, scalar_t, scalar_t>,   // normals ax, ay, az
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  bx, by, bz
                scalar_t, scalar_t, scalar_t>,   // normals bx, by, bz
            tuple<
                scalar_t, scalar_t, scalar_t,    // points  cx, cy, cz
                scalar_t, scalar_t, scalar_t> >  // normals cx, cy, cz
        const& tup) const // the function is const but data pointed to by pts
                          // is modified (thats the point of the function)
    {
        auto const& idxs = get<0>(tup);
        auto const& a_vals = get<1>(tup);
        auto const& b_vals = get<2>(tup);
        auto const& c_vals = get<3>(tup);

        if(get<0>(a_vals) != default_value)
        {
            int const a0 = get<0>(idxs) - offset;

            *(pts_x + a0) = get<0>(a_vals);
            *(pts_y + a0) = get<1>(a_vals);
            *(pts_z + a0) = get<2>(a_vals);
            *(nrs_x + a0) = get<3>(a_vals);
            *(nrs_y + a0) = get<4>(a_vals);
            *(nrs_z + a0) = get<5>(a_vals);
        }

        if(get<0>(b_vals) != default_value)
        {
            int const b0 = get<1>(idxs) - offset;

            *(pts_x + b0) = get<0>(b_vals);
            *(pts_y + b0) = get<1>(b_vals);
            *(pts_z + b0) = get<2>(b_vals);
            *(nrs_x + b0) = get<3>(b_vals);
            *(nrs_y + b0) = get<4>(b_vals);
            *(nrs_z + b0) = get<5>(b_vals);
        }

        if(get<0>(c_vals) != default_value)
        {
            int const c0 = get<2>(idxs) - offset;

            *(pts_x + c0) = get<0>(c_vals);
            *(pts_y + c0) = get<1>(c_vals);
            *(pts_z + c0) = get<2>(c_vals);
            *(nrs_x + c0) = get<3>(c_vals);
            *(nrs_y + c0) = get<4>(c_vals);
            *(nrs_z + c0) = get<5>(c_vals);
        }
    }
};

struct set_triangles
  : public unary_function<
        int,
        void>
{
    set_triangles(
        pointer<int> trs0, pointer<int> trs1, pointer<int> trs2,
        int const& nx, int const& ny, int const& nz,
        const_pointer<int> a0, const_pointer<int> b0, const_pointer<int> c0,
        const_pointer<uchar> cube_ids,
        const_pointer<int> tri_scan)
      : trs0(trs0), trs1(trs1), trs2(trs2),
        nx(nx), ny(ny), nz(nz),
        a0(a0), b0(b0), c0(c0),
        cube_ids(cube_ids),
        tri_scan(tri_scan)
    {}

    pointer<int> trs0, trs1, trs2;
    int const nx, ny, nz;
    const_pointer<int> a0, b0, c0;
    const_pointer<uchar> cube_ids;
    const_pointer<int> tri_scan;

    __platform__
    void
    operator()(int const& idx) const // const but modifies values pointed
                                     // to by trs0, trs1 and trs2
    {
        // only called for w = 0, ..., nw-2 where w is i, j or k.
        int const k = (idx / (nx-1)) / (ny-1);
        int const j = (idx / (nx-1)) % (ny-1);
        int const i = idx % (nx-1);

        uchar const& cube_id = *(cube_ids + idxer(i, j, k));

        if(cube_id == 0 || cube_id == 255)
        {
            return;
        }

        bool const* is_cut = mctable::isCut[cube_id];

        int global_idxs[12];

        if(is_cut[0])  global_idxs[0]  = *(a0 + idxer(i,   j,   k));
        if(is_cut[1])  global_idxs[1]  = *(b0 + idxer(i+1, j,   k));
        if(is_cut[2])  global_idxs[2]  = *(a0 + idxer(i,   j+1, k));
        if(is_cut[3])  global_idxs[3]  = *(b0 + idxer(i,   j,   k));
        if(is_cut[4])  global_idxs[4]  = *(a0 + idxer(i,   j,   k+1));
        if(is_cut[5])  global_idxs[5]  = *(b0 + idxer(i+1, j,   k+1));
        if(is_cut[6])  global_idxs[6]  = *(a0 + idxer(i,   j+1, k+1));
        if(is_cut[7])  global_idxs[7]  = *(b0 + idxer(i,   j,   k+1));
        if(is_cut[8])  global_idxs[8]  = *(c0 + idxer(i,   j,   k));
        if(is_cut[9])  global_idxs[9]  = *(c0 + idxer(i+1, j,   k));
        if(is_cut[10]) global_idxs[10] = *(c0 + idxer(i,   j+1, k));
        if(is_cut[11]) global_idxs[11] = *(c0 + idxer(i+1, j+1, k));


        int tri_idx = *(tri_scan + idxer(i, j, k));
        const char* tri_local_idx = mctable::caseTriangles[cube_id];
        for(int q = 0; tri_local_idx[q] != -1; q += 3)
        {
            *(trs0 + tri_idx) = global_idxs[tri_local_idx[q + 0]];
            *(trs1 + tri_idx) = global_idxs[tri_local_idx[q + 1]];
            *(trs2 + tri_idx) = global_idxs[tri_local_idx[q + 2]];

            tri_idx += 1;
        }
    }

    __platform__
    int
    idxer(int const& i, int const& j, int const& k) const
    {
        return k*nx*ny + j*nx + i;
    }
};

#endif
