/*
    struct Trim
    {
        Trim(t3<const int> n)
          : left(n.y*n.z, n.x),
            right(n.y*n.z, 0),
            edgeCases((n.x-1)*n.y*n.z),
            iter_helper(edgeCases.begin(), n.y*n.z, n.x-1)
        {}

        thrust::device_vector<int> left;
        thrust::device_vector<int> right;
        thrust::device_vector<uchar> edgeCases;
        iter_access_helper<uchar> iter_helper;

        using reference_type =
            thrust::tuple<
                int,
                int,
                typename thrust::device_vector<uchar>::iterator>;

        using iterator_tuple =
            thrust::tuple<
                typename thrust::device_vector<int>::iterator,
                typename thrust::device_vector<int>::iterator,
                typename iter_access_helper<uchar>::iterator>;

        using iterator = thrust::zip_iterator<iterator_tuple>;

        iterator begin()
        {
            return thrust::make_zip_iterator(
                thrust::make_tuple(
                    left.begin(),
                    right.begin(),
                    iter_helper.begin()));
        }

        iterator end()
        {
            return thrust::make_zip_iterator(
                thrust::make_tuple(
                    left.end(),
                    right.end(),
                    iter_helper.end()));
        }
    };
*/







/*
    struct set_trim_values
      : public thrust::binary_function<
            typename Trim::reference_type,                      // arg 1
            typename thrust::device_vector<scalar_t>::iterator, // arg 2
            typename Trim::reference_type>                      // out
    {
        set_trim_values(
            scalar_t const& isoval,
            int const& nx)
          : isoval(isoval),
            nx(nx)
        {}

        // This will not work on host because of using
        // edgeCases from device_vector...
        __host__
        typename Trim::reference_type
        operator()(
            typename Trim::reference_type trim_values,
            typename thrust::device_vector<scalar_t>::iterator curPoints)
        {
            // use curPoints to set xl, xr, curEdgeCases
            int xl = thrust::get<0>(trim_values);
            int xr = thrust::get<1>(trim_values);

            using iterator = typename thrust::device_vector<uchar>::iterator;
            iterator edgeCases = thrust::get<2>(trim_values);

            // TODO set all of isGE at once.

            bool isGE[2];
            isGE[0] = (curPoints[0] >= isoval);
            for(int i = 1; i != nx; ++i)
            {
                isGE[i%2] = (curPoints[i] >= isoval);

//                edgeCases[i-1] = calcCaseEdge(isGE[(i+1)%2], isGE[i%2]);
//
                if(*(edgeCases + i-1) == 1 || *(edgeCases + i-1) == 2)
                {
                    if(xl > xr)
                    {
                        xl = i-1;
                    }
                    xr = i;
                }
            }

            return thrust::make_tuple(xl, xr, edgeCases);
        }

        scalar_t const& isoval;
        int const& nx;
    };

    void pass1()
    {
        thrust::transform(
            trim.begin(),       // input1
            trim.end(),
            image.ray_begin(),  // input2
            trim.begin(),       // output
            set_trim_values(    // binary function (*input1, *input2)
                isoval,
                n.x));
    }
*/

