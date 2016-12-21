//
// Part of the observer/subject pattern
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "subject.hh"
#include "observer.hh"
#include "stream_manager.hh"

void 
fluke::Subject::attach( Observer *o ) {
    obs_.push_back( o );
}

void 
fluke::Subject::detach( Observer *o ) {
    // observer 'o' exists in obs_
    obs_.erase( std::find( obs_.begin(), obs_.end(), o ) );
}

void 
fluke::Subject::detachAll() {
    obs_.clear();
}

void 
fluke::Subject::notify() 
{ notify( this ); }

void
fluke::Subject::notify( Subject *s ) {
    for( obs_iter i = obs_.begin(); i != obs_.end(); ++i ) {
        ( **i ).update( s );
    }
}

