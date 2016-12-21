//
// Agent identification tag
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_AGENT_TAG_
#define _FLUKE_AGENT_TAG_

#include "defs.hh"

namespace fluke {
    
    /// \class AgentTag
    /// \brief Unique identification for agents
    /// 
    /// Ancestor tracing of agents is done with a small \c class named
    /// \c AgentTag.
    class AgentTag {
        public:
        /// Constructor for convenience
        AgentTag() : time( -1 ), x( -1 ), y( -1 ), i( -1 ) {}
        AgentTag( long t, int x, int y, int i ) 
            : time( t ), x( x ), y( y ), i( i ) {}
        ~AgentTag() {}
        
        /// As a string
        std::string str() const;

        public:        
        /// Time of birth
        long time;
        /// x coordinate of its birth location
        int x;
        /// y coordinate of its birth location
        int y;
        /// index for those agents that are born in the same timestep as their
        /// parent
        int i;
    };
    
    inline bool operator==( const AgentTag &t1, const AgentTag &t2 ) 
    { return t1.time == t2.time && t1.x == t2.x 
        && t1.y == t2.y && t1.i == t2.i; }
    
    inline std::string AgentTag::str() const
    { return boost::lexical_cast< std::string >( time ) + "-" +
        boost::lexical_cast< std::string >( x ) + "-" +
        boost::lexical_cast< std::string >( y ) + "-" +
        boost::lexical_cast< std::string >( i ); }

}
#endif
