//
// Implementation of a simple agent (not complete).
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "simple_agent.hh"
#include "population.hh"


float fluke::SimpleAgent::birth_rate_ = 0.0;
float fluke::SimpleAgent::death_rate_ = 0.0;


fluke::SimpleAgent::SimpleAgent() : Agent(), score_( birth_rate_ ) {}

fluke::SimpleAgent::SimpleAgent( const SimpleAgent &ag ) : Agent() {
    SimpleAgent::copy( ag );
}

fluke::Agent* 
fluke::SimpleAgent::clone() const {
    return new SimpleAgent( *this );
}

void 
fluke::SimpleAgent::copy( const Agent &ag ) {
    const SimpleAgent &sag = dynamic_cast< const SimpleAgent & >( ag );
    Agent::copy( ag );
    score_ = sag.score_;
}

void 
fluke::SimpleAgent::step( Population &pop ) {
    float rr = uniform();
    if( rr < death_rate_ ) {
        dying_ = true;
    }
}

fluke::Agent*
fluke::SimpleAgent::sibling() {
    Agent* result;
    result = this->clone();
    //SimpleAgent *aux = dynamic_cast< SimpleAgent* >( result );
    return result;
}

void 
fluke::SimpleAgent::write( std::ostream &os ) const {
    os << "<agent birth=\"" << me_.time << "\" x=\"" << me_.x;
    os << "\" y=\"" << me_.y << "\" i=\"" << me_.i << "\"";
    os << ">\n<class>SimpleAgent</class>\n</agent>\n";
}

void
fluke::SimpleAgent::birthRate( float f )
{ birth_rate_ = f; }

float
fluke::SimpleAgent::birthRate()
{ return birth_rate_; }

void
fluke::SimpleAgent::deathRate( float f )
{ death_rate_ = f; }

float
fluke::SimpleAgent::deathRate()
{ return death_rate_; }

