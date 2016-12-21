//
// Adaptor for passing two agents to an observer
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_DUO_AGENT_
#define _FLUKE_DUO_AGENT_

#include "defs.hh"

namespace fluke {
 
    struct DuoAgent : public Subject, std::pair< Agent*, Agent * > {
        DuoAgent() : Subject(), std::pair< Agent*, Agent* >() {};
        DuoAgent( Agent *a, Agent *b ) 
            : Subject(), std::pair< Agent*, Agent* >( a, b ) {};
        ~DuoAgent() {};
    };
    
    inline DuoAgent make_duo( Agent *a, Agent *b ) 
    { return DuoAgent( a, b ); }
}
#endif
