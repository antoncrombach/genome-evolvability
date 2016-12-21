//
// XML SAX reader for agents
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _AGENT_READER_
#define _AGENT_READER_

#include "defs.hh"
#include "bsite.hh"
#include "transfac.hh"
#include "module_dstream.hh"
#include "ordinary_dstream.hh"
#include "retroposon.hh"
#include "repeat.hh"
#include "chromosome.hh"
#include "genome.hh"
#include "simple_agent.hh"
#include "module_agent.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace fluke {

// Xerces shit.
XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    /// \class AgentReader
    /// \brief Read in an agent from an xml file.
    ///
    /// AgentReader reads a genome from an xml file. If multiple genomes are
    /// present (I think) it reads just the first agent and returns it.
    class AgentReader : public DefaultHandler {
        public:
            /// Constructor.
            AgentReader( Factory * );
            /// Destructor.
            ~AgentReader() {};
            
            /// Read opening xml tag.
            void startElement( const XMLCh* const uri, 
                const XMLCh* const localname, const XMLCh* const qname,
                const Attributes &attrs );
            /// Read tag contents. 
            void characters( const XMLCh* const chars, const uint length );
            /// Read closing xml tag.
            void endElement( const XMLCh* const uri, 
                const XMLCh* const localname, const XMLCh* const qname );
            /// 
            void resetDocument();
        
            // A few functions for catching warnings and errors
            void warning( const SAXParseException &exc );
            void error( const SAXParseException &exc );
            void fatalError( const SAXParseException &exc );
            void resetErrors();
    
            bool sawErrors() const;
            bool done() const;
            
            // Return the agent
            Agent* agent();
            
            // Set the type
            void type( int );
        
        private:
            Factory *factory_;
            bool sawErrors_, done_, class_, ess_, mods_;
            
            std::string agent_class_;
            int type_;
            
            std::list< ChromosomeElement* > *chr_;
            Chromosome *chromo_;
            Genome *genome_;
            Agent *agent_;
            AgentTag me_, mother_;
    };

inline bool AgentReader::sawErrors() const
{ return sawErrors_; }

inline bool AgentReader::done() const
{ return done_; }

inline void AgentReader::type( int tt )
{ type_ = tt; }

}
#endif
