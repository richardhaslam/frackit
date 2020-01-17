#include <frackit/geometry/segment.hh>
#include <frackit/magnitude/length.hh>
#include <frackit/intersection/intersect.hh>

enum IntersectionType { point, segment, empty };

template<class G>
void checkResultGeometry(const G& geometry, IntersectionType expected)
{
    throw std::runtime_error(std::string("Unexpected intersection geometry"));
}

template<int wd>
void checkResultGeometry(const Frackit::EmptyIntersection<wd>& empty, IntersectionType expected)
{
    std::cout << "Found empty intersection" << std::endl;
    if (expected != IntersectionType::empty)
        throw std::runtime_error(std::string("Unexpected empty intersection"));
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Point<CT, wd>& p, IntersectionType expected)
{
    std::cout << "Found intersection point at " << p << std::endl;
    if (expected != IntersectionType::point)
        throw std::runtime_error(std::string("Got an unexpected point intersection"));
    std::cout << "Test passed" << std::endl;
}

template<class CT, int wd>
void checkResultGeometry(const Frackit::Segment<CT, wd>& segment, IntersectionType expected)
{
    std::cout << "Found intersection segment with corners "
              << segment.source() << " - " << segment.target()
              << " and length " << computeLength(segment) << std::endl;
    if (expected != IntersectionType::segment)
        throw std::runtime_error(std::string("Got an unexpected segment intersection"));
    std::cout << "Test passed" << std::endl;
}

//! test segment-segment intersections
int main()
{
    using ctype = double;
    using Segment3 = Frackit::Segment<ctype, 3>;
    using Point3 = typename Segment3::Point;

    std::vector<ctype> scales({1e-5, 1, 1e5});

    for (auto f : scales)
    {
        std::cout << "Checking scale factor " << f << std::endl;

        // check segments in x-y plane
        Segment3 s1( Point3(0.0,   0.0,   0.0), Point3(1.0*f, 1.0*f, 0.0) );
        Segment3 s2( Point3(0.5*f, 0.5*f, 0.0), Point3(1.0*f, 1.0*f, 0.0) );
        Segment3 s3( Point3(0.0,   1.0*f, 0.0), Point3(1.0*f, 1.0*f, 0.0) );
        Segment3 s4( Point3(0.0,   1.0*f, 0.0), Point3(1.0*f, 0.0,   0.0) );

        auto result = intersect(s1, s2);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);

        result = intersect(s1, s3);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        result = intersect(s1, s4);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        // check segments in x-y plane
        Segment3 s5( Point3(0.0,   0.0,   0.0),   Point3(1.0*f, 1.0*f, 1.0*f) );
        Segment3 s6( Point3(0.5*f, 0.5*f, 0.5*f), Point3(1.0*f, 1.0*f, 1.0*f) );
        Segment3 s7( Point3(0.0,   1.0*f, 0.0),   Point3(1.0*f, 1.0*f, 1.0*f) );
        Segment3 s8( Point3(0.0,   1.0*f, 0.0),   Point3(1.0*f, 0.0,   1.0*f) );

        result = intersect(s5, s6);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::segment); }, result);

        result = intersect(s5, s7);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        result = intersect(s5, s8);
        std::visit([&] (auto&& is) { checkResultGeometry(is, IntersectionType::point); }, result);

        std::cout << "All tests passed"  << std::endl;
    }

    return 0;
}
