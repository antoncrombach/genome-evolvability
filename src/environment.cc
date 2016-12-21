//
// Implementation of three environments.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "environment.hh"

fluke::Environment::Environment() 
    : model_( 0 ), generator_env_( 18 ), uniform_env_( generator_env_ ) {}

void
fluke::Environment::initialise( uint s ) {
    generator_env_.seed( s );
    uniform_gen_type aux( generator_env_ );
    uniform_env_ = aux;
    
    if( model_ != 0 ) {
        model_->population().evaluate( *this );
    }
}

//
// Constant environment
//
fluke::ConstantEnvironment::ConstantEnvironment( int k ) 
    : Environment(), copies_( k, 0 ) {}

//
// Periodic environment
//
fluke::PeriodicEnvironment::PeriodicEnvironment( int k ) 
    : Environment(), copies_( k, 0 ), periods_( k, 0 ), offsets_( k, 0 ),
      state_a_( k, 0 ), state_b_( k, 0 ) {}

void
fluke::PeriodicEnvironment::states( int k, int low, int high ) {
    state_a_[ k ] = low;
    state_b_[ k ] = high;
}

void
fluke::PeriodicEnvironment::initialise( uint s ) {
    // ignoring s
    for( uint i = 0; i != copies_.size(); ++i ) {
        copies_[ i ] = state_a_[ i ];
    }
}

void
fluke::PeriodicEnvironment::fluctuate( long time ) {
    bool aux = false;
    for( uint i = 0; i != copies_.size(); ++i ) {
        if( ( time - offsets_[ i ] ) % periods_[ i ] == 0 ) {
            aux = true;
            copies_[ i ] = ( copies_[ i ] != state_a_[ i ] ) ?
                state_a_[ i ] : state_b_[ i ];
            // inform observers
            notify();
        }
    }
    if( aux ) {
        model_->population().evaluate( *this );
    }
}

//
// Poisson environment
//
fluke::PoissonEnvironment::PoissonEnvironment( int k ) 
    : Environment(), copies_( k, 0 ), lambdas_( k, 0.0 ),
      state_a_( k, 0 ), state_b_( k, 0 ) {}

void
fluke::PoissonEnvironment::states( int k, int low, int high ) {
    state_a_[ k ] = low;
    state_b_[ k ] = high;
}

void
fluke::PoissonEnvironment::initialise( uint seed ) {
    for( uint i = 0; i != copies_.size(); ++i ) {
        copies_[ i ] = state_a_[ i ];
    }
    Environment::initialise( seed );
}

void
fluke::PoissonEnvironment::fluctuate( long time ) {
    bool aux = false;
    for( uint i = 0; i != copies_.size(); ++i ) {
        if( uniform_env_() < lambdas_[ i ] ) {
            aux = true;
            copies_[ i ] = ( copies_[ i ] != state_a_[ i ] ) ?
                state_a_[ i ] : state_b_[ i ];
            // inform observers
            notify();
        }
    }
    if( aux ) {
        model_->population().evaluate( *this );
    }
}

