//
// Environment base class and childs with simple (periodic) changes.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_ENVIRONMENT_H_
#define _FLUKE_ENVIRONMENT_H_

#include "defs.hh"
#include "subject.hh"
#include "population.hh"

namespace fluke {

    /// \class Environment
    /// \brief Abstract base class for modelling the environment.
    ///
    /// The interplay between the environment and the population may become
    /// quite intricate. The current scope of the environment is that it is
    /// some external mechanism influencing the population. The population has
    /// no influence on the environment though.
    class Environment : public Subject {
        public:
        /// Destructor (empty).
        virtual ~Environment() {};
        
        /// Initialise the environment
        virtual void initialise( uint );
        /// The environment changes, the population is alerted
        virtual void fluctuate( long ) = 0;
        /// Get per module the number of `optimal' gene copies
        virtual int expectedCopies( int ) const = 0;
        
        /// Set a reference to the model
        void model( Model * );
        /// Get a reference to the model
        Model* model() const;
        
        protected:
        /// Hidden constructor
        Environment();
        
        protected:
        /// Reference to the model
        Model *model_;
        /// The environment has its own uniform random nr generator
        base_generator_type generator_env_;
        uniform_gen_type uniform_env_;
    };
    
    inline void Environment::model( Model *m )
    { model_ = m; }
    
    inline Model* Environment::model() const
    { return model_; }
    
    
    /// \class ConstantEnvironment
    /// \brief Most simple environment, the not-changing one.
    ///
    /// The constant environment is used as a reference to see the influence 
    /// of the other environments.
    class ConstantEnvironment : public Environment {
        public:
        /// Constructor
        ConstantEnvironment( int );
        /// Destructor
        ~ConstantEnvironment() {}
        
        /// No fluctuation in a constant environment.
        virtual void fluctuate( long ) {}
        /// What is optimal...
        virtual int expectedCopies( int ) const;
        /// Give for a module, the nr of copies of the genes needed to be 
        /// well adapted to this environment
        void expectedCopies( int, int );
        
        private:
        std::vector< int > copies_;
    };

    inline int ConstantEnvironment::expectedCopies( int k ) const
    { return copies_[ k ]; }

    inline void ConstantEnvironment::expectedCopies( int k, int c )
    { copies_[ k ] = c; }

    /// \class PeriodicEnvironment
    /// \brief Environment changes periodically
    ///
    /// Now we know when the environment changes. (Contrast to 
    /// PoissonEnvironment.)
    class PeriodicEnvironment : public Environment {
        public:
        /// Constructor
        PeriodicEnvironment( int );
        /// Destructor
        ~PeriodicEnvironment() {}
        
        /// Initialise every module to the \c low number of gene copies
        virtual void initialise( uint );
        /// Each module had a chance of switching the desired number of
        /// gene copies. The change is a toggle between two states (values)
        void fluctuate( long );
        /// What is optimal...
        virtual int expectedCopies( int ) const;
        /// Give for a module, the nr of copies of the genes needed to be 
        /// well adapted to this environment
        void expectedCopies( int, int );
        /// Set the length of a period.
        void period( int, long );
        /// Set the two allowed states between which the environment 
        /// changes a module. Usually the two states are refered to as 
        /// \c low and \c high.
        void states( int, int, int );
        /// Periods may be shifted to overlap more or less.
        void offset( int, long );
        
        private:
        std::vector< int > copies_;
        std::vector< long > periods_;
        std::vector< long > offsets_;
        
        std::vector< int > state_a_;
        std::vector< int > state_b_;
    };

    inline int PeriodicEnvironment::expectedCopies( int k ) const 
    { return copies_[ k ]; }

    inline void PeriodicEnvironment::period( int k, long l ) 
    { periods_[ k ] = l; }

    inline void PeriodicEnvironment::offset( int k, long l ) 
    { offsets_[ k ] = l; }

    /// \class PoissonEnvironment
    /// \brief Environment changes according to stochastic Poisson process.
    ///
    /// As a short cut to expression levels we take number of gene copies as
    /// the criterium for being adapted to the environment. 
    class PoissonEnvironment : public Environment {
        public:
        /// Constructor
        PoissonEnvironment( int );
        /// Destructor
        ~PoissonEnvironment() {}
        
        /// Initialise every module to the \c low number of gene copies
        virtual void initialise( uint );
        /// Each module had a chance of switching the desired number of
        /// gene copies. The change is a toggle between two states (values)
        void fluctuate( long );
        /// Get the current `optimal' number of gene copies
        virtual int expectedCopies( int ) const;
        /// Set the probability of toggling between two states
        void lambda( int, double );
        /// Set the two allowed states between which the environment 
        /// changes a module. Usually the two states are refered to as 
        /// \c low and \c high.
        void states( int, int, int );
        
        private:
        std::vector< int > copies_;
        std::vector< double > lambdas_;
        
        std::vector< int > state_a_;
        std::vector< int > state_b_;
    };
    
    inline int PoissonEnvironment::expectedCopies( int k ) const 
    { return copies_[ k ]; }

    inline void PoissonEnvironment::lambda( int k, double l ) 
    { lambdas_[ k ] = l; }

}
#endif
