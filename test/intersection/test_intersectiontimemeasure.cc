#include <cmath>
#include <chrono>

#include <frackit/geometry/disk.hh>
#include <frackit/intersection/intersect.hh>
#include <frackit/occ/breputilities.hh>
#include <frackit/common/math.hh>

//! test that stops the average time needed for intersections
//! can be used together with git bisect to check if certain
//! changes affected the efficiency of the intersection computations
int main()
{
    using namespace Frackit;
    using ctype = double;
    using Disk = Disk<ctype>;
    using Point = typename Disk::Point;
    using Direction = typename Disk::Direction;
    using Vector = typename Direction::Vector;

    std::vector<ctype> scales({1.0e-5, 1, 1e5});
    for (auto f : scales)
    {
        std::cout << "Scale factor " << f << std::endl;

        const Direction e1(Vector(1.0, 0.0, 0.0));
        const Direction e2(Vector(0.0, 1.0, 0.0));
        const Disk disk(Point(0.0, 0.0, 0.0), e1, e2, f, f);
        const auto diskShape = OCCUtilities::getShape(disk);

        const Disk disk2(Point(0.0, 0.0, 10.0*f), e1, e2, f, f);
        const auto diskShape2 = OCCUtilities::getShape(disk2);

        using std::sin;
        using std::cos;
        using std::abs;

        // create 1000 disks at different angles
        static constexpr std::size_t numTests = 1000;
        std::vector<Disk> testDisks; testDisks.reserve(numTests);
        std::vector<TopoDS_Face> testDiskShapes; testDiskShapes.reserve(numTests);

        static constexpr ctype deltaAngle = M_PI/(numTests+1);
        for (std::size_t i = 1; i <= numTests; ++i)
        {
            ctype angle = i*deltaAngle;
            const Direction majAxis( Vector(cos(angle), 0.0, sin(angle)) );

            testDisks.emplace_back(Point(0.0, 0.0, 0.0), majAxis, e2, 0.5*f, 0.5*f);
            testDiskShapes.emplace_back(OCCUtilities::getShape(testDisks.back()));
        }

        // test non-empty intersections
        std::cout << std::endl;
        auto start = std::chrono::steady_clock::now();
        for (const auto& d : testDisks)
            intersect(disk, d);
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration<double>(end-start).count();
        std::cout << "Non-empty intersections with internal geometry took "
                  << duration << " -> " << duration/numTests << " on average" << std::endl;

        start = std::chrono::steady_clock::now();
        for (const auto& d : testDiskShapes)
            intersect(diskShape, d);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end-start).count();
        std::cout << "Non-empty intersections with shape representation took "
                  << duration << " -> " << duration/numTests << " on average" << std::endl;

        // test empty intersections
        std::cout << std::endl;
        start = std::chrono::steady_clock::now();
        for (const auto& d : testDisks)
            intersect(disk2, d);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end-start).count();
        std::cout << "Empty intersections with internal geometry took "
                  << duration << " -> " << duration/numTests << " on average" << std::endl;

        start = std::chrono::steady_clock::now();
        for (const auto& d : testDiskShapes)
            intersect(diskShape2, d);
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end-start).count();
        std::cout << "Empty intersections with shape representation took "
                  << duration << " -> " << duration/numTests << " on average" << std::endl;
    }

    return 0;
}
