//
// Interface of an abstract agent.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_AGENT_H_
#define _FLUKE_AGENT_H_

#include "defs.hh"
#include "subject.hh"
#include "agent_tag.hh"

namespace fluke {

    /// \class Agent 
    /// \brief Abstract agent. 
    /// 
    /// The minimal interface of any \c Agent is defined. It provides methods 
    /// for cloning agents, updating, reproducing and it signals when it is 
    /// dying.
    /// In addition a standard way of writing into an \c ostream is present.
    class Agent : public Subject {
        public:
        /// Virtual destructor.
        virtual ~Agent() {};
        /// Cloning an agent. A flexible and easy way of producing copies
        /// of an agent. Also known as a prototyping OOP pattern.
        virtual Agent* clone() const = 0;
        /// Copy an agent.
        virtual void copy( const Agent & );
        
        /// Initialise the agent.
        virtual void initialise() = 0;
        /// Perform all actions of a timestep. The method needs to be
        /// overriden in child classes to obtain more interesting behaviour.
        virtual void step( Population & ) = 0;
        /// Create a sibling. The \c fluke equivalent of reproduction.
        virtual Agent* sibling() = 0;
        
        /// Signature for score get-method.
        virtual double score() const = 0;
        /// Signature for calculation of the score.
        virtual void evaluate( const Environment & ) = 0;
        
        /// Signature for writing to output streams.
        virtual void write( std::ostream & ) const = 0;
        /// Is the agent dying?
        bool dying() const;
        
        /// Set the ID of the immediate ancestor.
        void parentTag( AgentTag );            
        /// Get the ID of the direct ancestor.
        AgentTag parentTag() const;
        /// Set the agent's ID.
        void myTag( AgentTag );
        /// Get the agent's ID.
        AgentTag myTag() const;
        /// Set type of agent (used in competition experiments).
        void type( int );
        /// Get type of agent.
        int type() const;
        
        /// Get its (genotypical) distance.
        virtual int distance() const = 0;
            
        protected:
        /// Hidden constructor.
        Agent();
        /// Hidden constructor with agent ID-tag.
        Agent( AgentTag );
        /// Hidden constructor with agent type.
        Agent( int );
        /// Hidden copy constructor.
        Agent( const Agent & );
            
        protected:
        /// Ready to be erased from the population (flag).
        bool dying_;
        /// My ID-tag.
        AgentTag me_;
        /// Parent's ID-tag.
        AgentTag ancestor_;
        /// Type (competition experiments).
        int type_;
    };

    /// Overloaded \c << operator for easy writing to streams (xml)
    inline std::ostream& operator<<( std::ostream& os,
        const Agent& ag )
    { ag.write( os ); return os; }

    inline bool Agent::dying() const
    { return dying_; }

    inline void Agent::parentTag( AgentTag t )
    { ancestor_ = t; }
    
    inline AgentTag Agent::parentTag() const
    { return ancestor_; }

    inline void Agent::myTag( AgentTag t )
    { me_ = t; }
    
    inline AgentTag Agent::myTag() const
    { return me_; }
    
    inline void Agent::type( int t )
    { type_ = t; }
    
    inline int Agent::type() const
    { return type_; }
}
#endif

