//
// Interface of an simple agent.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_SIMPLEAGENT_H_
#define _FLUKE_SIMPLEAGENT_H_

#include "defs.hh"
#include "agent.hh"

namespace fluke {

    /// \class SimpleAgent
    /// \brief Dummy agent that is used for testing \c Population.
    ///
    /// Spatial embedding and population dynamics are an important part of 
    /// the model. Most GAs and all ODEs ignore space, yet it has been shown 
    /// multiple times that the influence of space on the behaviour of the
    /// model is not trivial.
    class SimpleAgent : public Agent {
        public:
            /// Constructor.
            SimpleAgent();
            /// Copy constructor.
            SimpleAgent( const SimpleAgent & );
            /// Destructor.
            virtual ~SimpleAgent() {}
            /// Cloning an agent.
            virtual Agent* clone() const;
            /// Copy the agent into \c this
            virtual void copy( const Agent & );

            /// Dummy (FIX)
            virtual void initialise() {};
            /// Perform one timestep. It consists of checking whether the
            /// agent dies or survives another timestep.
            virtual void step( Population & );
            /// Reproduce. The sibling has its age reset to zero.
            virtual Agent* sibling();

            /// Get the score of a \c SimpleAgent, equal to its birth rate.
            virtual double score() const;
            /// Empty function
            virtual void evaluate( const Environment & );
            /// Another empty function
            virtual int distance() const;            
            /// Write a text representation of the agent to an output stream.
            virtual void write( std::ostream & ) const;
            
        public:
            /// Set the birth rate.
            static void birthRate( float );
            /// Get the birth rate.
            static float birthRate();
            /// Set the death rate.
            static void deathRate( float );
            /// Get the death rate.
            static float deathRate();
            
        protected:
            /// Score (birth rate) of a \c SimpleAgent.
            float score_;

        protected:
            /// Birth rate.
            static float birth_rate_;
            /// Death rate.
            static float death_rate_;
    };

    inline double SimpleAgent::score() const
    { return score_; }

    inline void SimpleAgent::evaluate( const Environment &env )
    { throw "Not implemented for SimpleAgent"; }
    
    inline int SimpleAgent::distance() const
    { throw "Not implemented for SimpleAgent"; }
}
#endif

