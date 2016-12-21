//
// Definitions of globals needed by entire project.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_DEFS_H_
#define _FLUKE_DEFS_H_

#include <unistd.h>

#include <map>
#include <list>
#include <cmath>
#include <string>
#include <stack>
#include <vector>
#include <cassert>
#include <limits>
#include <utility>
#include <numeric>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <ext/numeric>

#include <boost/regex.hpp>
#include <boost/utility.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/lagged_fibonacci.hpp>

//#include <boost/vector_property_map.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graphviz.hpp>

#include "util.hh"
#include "pool.hh"

//
// Memory policy for arguments and return values: 
// - expect aliases (pointers and references) as arguments. Methods need to 
//   copy things themselves.
// - expect aliases as return values.
//
// Option: do not return objects, instead pass them as a reference argument.
// Removes the hazard of implicitly copying objects.

#ifdef DEBUG
using std::cout;
using std::endl;
#endif

// Random number generators
typedef boost::lagged_fibonacci607 base_generator_type;
typedef boost::variate_generator< base_generator_type, boost::uniform_int<> > 
    randrange_gen_type;
typedef boost::uniform_01< base_generator_type > uniform_gen_type;

/// \namespace fluke The project is placed in the namespace \c fluke.
/// A \em fluke is a stroke of good luck. At the start (July 2004) the whole 
/// project seemed quite adventurous to me, hence the name. Other namespaces 
/// used (only explicitly though) are \c std and \c boost.
///
/// Some design decisions that have been made (and may be reconsidered) are 
/// listed below.
/// - rethink CA vs agent-based modelling: for the moment we keep CA-like 
///   modelling
/// - for the moment ignoring binding sites (gene duplication may interrupt
///   other upstream regions, in contrast to Otto)
/// - the id of an agent is (time of birth, x, y, i), where (x, y) is its grid
///   location and i is the index for multiple divisions in 1 timestep
/// - mutations are maybe not completely in a biological order, but we can
///   always change that quite easily
/// - distinction between lifetime and replication mutations is not in place
/// - using xerces-c SAX2 parser for xml access of XML agent and population
///   reader
/// - xml is only nice for saving a simulation or whole genomes 
/// - extra options needed for g++ to be able to write >2Gb files
///   (-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE64_SOURCE)    
///   (fixed in libstd6.0+)
/// - new fitness function: exp( -x / 23.0 )
/// - new score function: distance to total nr of wanted genes
namespace fluke {

    /// Version number of the application in the format "a.b.c". \c a is 
    /// increased when either major changes are made to the program or a fully
    /// stable, 'public' release is made. \c b signals the inclusion of new
    /// features in the current release. \c c is used for series of bugfixes.
    const std::string VERSION = "2.8.1";

    /// Labels of short sequences are integers.
    typedef int label;
    /// unsigned int is tooo much typing
    typedef unsigned int uint;

    /// Coz ya always need the little fella
    const double PI = 3.14159265;

    // all the classes declared
    class ShortSeqManager;
    class ChromosomeElement;
    class BindingSite;
    class Downstream;
    class TranscriptionFactor;
    class ModuleDownstream;
    class OrdinaryDownstream;
    class Retroposon;
    class Repeat;
    class Chromosome;
    class Genome;
    class MutateRates;
    // reviewer 2
    class Centromere;

    class Agent;    
    class AgentTag;
    class DuoAgent;
    class SimpleAgent;
    class ModuleAgent;

    class Environment;
    class ConstantEnvironment;
    class PoissonEnvironment;
    class PeriodicEnvironment;
    
    struct Location;
    class Population;
    class ScalingScheme;
    class NoScaling;
    class LinearScaling;
    class SelectionScheme;
    class RandSelection;
    class ProbalisticSelection;

    // rethink with concept of owner and being observed uncoupled
    class Observer;
    class LogObserver;
    class AsyncLogObserver;
    class Subject;

    class LogCsvMutations;
    class LogCsvGenes;
    class LogCsvPrunedDist;
    class LogCsvDistances;
    class LogCsvScores;
    class LogCsvEnvironment;
    class LogCsvAncestors;
    class LogXmlGenomes;
    class LogXmlAgentTrace;
    class LogXmlEnvGenomes;
    class LogPopulationSize;
    class LogCsvPopulationDistances;
    
    class Factory;
    class AgentReader;
    class PopulationReader;
    class ObserverManager;
    class StreamManager;
    class Model;
    class Config;
    class Fluke;

    /// Random number generator. It is a global object as it is used 
    /// throughout the entire program. It needs to be used in conjunction
    /// with a distribution.
    extern base_generator_type generator;
    /// Uniform random numbers [0,1). It is a global object to provide easy
    /// access in the entire program.
    extern uniform_gen_type uniform;

    /// Generate random number of type T from the interval \f$[0,n)\f$.
    template< class T > T
    rand_range( T n ) {
        return static_cast< T >( n * uniform() );
    }
}
#endif

