
#if !defined(_LIB_QDT_WRAPPER_HXX_)
#define _LIB_QDT_WRAPPER_HXX_

#include <array>
#include <vector>

namespace wrapper
{
    struct Triangulation
    {
        std::vector<std::array<size_t, 3>> triangles;
        std::vector<size_t> rep;
        int unique_points = 0;
    };
}

#endif
