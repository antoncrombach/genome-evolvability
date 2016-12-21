//
// Implementation of abstract agent.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "agent.hh"

fluke::Agent::Agent() 
: dying_( false ), me_(), ancestor_(), type_( -1 ) {}

fluke::Agent::Agent( AgentTag t ) 
: dying_( false ), me_( t ), ancestor_(), type_( -1 ) {}

fluke::Agent::Agent( int tt ) 
: dying_( false ), me_(), ancestor_(), type_( tt ) {}

fluke::Agent::Agent( const Agent &ag ) {
    copy( ag );
}

void
fluke::Agent::copy( const Agent &ag ) {
    dying_ = ag.dying_;
    me_ = ag.me_;
    ancestor_ = ag.ancestor_;
    type_ = ag.type_;
}

