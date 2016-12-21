//
// Mutating parameters can be done in different ways.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_MUTATE_RATES_
#define _FLUKE_MUTATE_RATES_

#include "defs.hh"

namespace fluke {
    class MutateRates {
        public:
        virtual ~MutateRates() {};
        virtual void mutate( Chromosome & ) = 0;
        
        protected:
        MutateRates() {};
    };
    
    class FixedMutateRates : public MutateRates {
        public:
        FixedMutateRates() : MutateRates() {};
        virtual ~FixedMutateRates() {};
        virtual void mutate( Chromosome & );
    };
    
    class LinearMutateRates : public MutateRates {
        public:
        LinearMutateRates() : MutateRates() {};
        virtual ~LinearMutateRates() {};
        virtual void mutate( Chromosome & );
    };
    
    class UniformMutateRates : public MutateRates {
        public:
        UniformMutateRates() : MutateRates(), high_( 1.0 ), low_( 0.0 ) {};
        UniformMutateRates( double ll, double hh ) 
            : MutateRates(), high_( hh ), low_( ll ) {};
        virtual ~UniformMutateRates() {};
        virtual void mutate( Chromosome & );
        
        private:
        double high_, low_;
    };
}
#endif
