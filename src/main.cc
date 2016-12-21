//
// Starting point of the fluke program.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "defs.hh"
#include "fluke.hh"

#ifdef DEBUG
#include <boost/timer.hpp>
#endif

using namespace std;
using namespace fluke;

// globals...
base_generator_type fluke::generator( 18 );
uniform_gen_type fluke::uniform( fluke::generator );

int
main( int argc, char **argv ) {
    int result = 1;
    try {
#ifdef DEBUG
        boost::timer tt;
#endif
        Fluke f( argc, argv );
        result = f.run();
#ifdef DEBUG
        std::cout << "Timing: " << tt.elapsed() << " seconds\n";
#endif
    } catch( char *e ) {
        std::cout << "Exception: " << e << std::endl;
    } catch( exception &e ) {
        std::cout << "Exception: " << e.what() << std::endl;
    }
    return result;
}

