// apps/probes/hello/main.cpp -- build/link/license smoke test (checks no geometry correctness).
// Verifies public headers resolve, sggk + ModelingPrim link, and a cylinder can be created.
// Exit codes: 0 pass, 1 body creation failed (license/link), 2 unexpected exception.

#include <Foundation/init.h>
#include <Foundation/SGGK_Version.h>
#include <ModelingPrim/API.h>
#include <iostream>

using namespace sggk;

// Smoke entry: init -> print version -> make cylinder -> fini; result encoded in the exit code.
int main()
{
    try
    {
        init();

        // CI logs the version string to confirm the same build artifact.
        std::cout << "[VERSION] " << SGGK_VERSION_COMPLETE
                  << "  ext=" << SGGK_VERSION_STRING_EXT << std::endl;

        // Cylinder creation also verifies the license and ModelingPrim linkage:
        // default Ucs3D, radius 200 mm, height 500 mm, full circle PI2 (2*pi rad).
        auto ret = api_make_solid_cylinder(Ucs3D(), 200.0, 500.0, PI2, true);
        if (ret == nullptr || !ret->Succeeded())
        {
            std::cerr << "[FAIL] api_make_solid_cylinder" << std::endl;
            // fini() is called on every path, including this failure return.
            fini();
            return 1;
        }
        std::cout << "[SMOKE] solid cylinder OK" << std::endl;

        fini();
    }
    catch (...)
    {
        std::cerr << "[FAIL] unhandled exception" << std::endl;
        return 2;
    }
    return 0;
}
