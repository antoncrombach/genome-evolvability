//
// XML SAX reader for populations
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _POPULATION_READER_
#define _POPULATION_READER_

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
#include "population.hh"

#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/Attributes.hpp>

namespace fluke {

XERCES_CPP_NAMESPACE_USE

XERCES_CPP_NAMESPACE_BEGIN
    class AttributeList;
XERCES_CPP_NAMESPACE_END

    class PopulationReader : public DefaultHandler {
        public:
            /// Constructor
            PopulationReader( Factory * );
            /// Destructor
            ~PopulationReader() {};
            
            void startElement( const XMLCh* const uri, 
                const XMLCh* const localname, const XMLCh* const qname,
                const Attributes &attrs );
            void characters( const XMLCh* const chars, 
                const unsigned int length );
            void endElement( const XMLCh* const uri, 
                const XMLCh* const localname, const XMLCh* const qname );     
            void resetDocument();
        
            // A few functions for catching warnings and errors
            void warning( const SAXParseException &exc );
            void error( const SAXParseException &exc );
            void fatalError( const SAXParseException &exc );
            void resetErrors();
    
            bool sawErrors() const;
            bool done() const;
            void type( int );
            
            // Return the population
            boost::tuple< std::vector< Agent* >, std::vector< Location > >
                population();
        
        private:
            Factory *factory_;
            bool sawErrors_, done_, class_;
            
            std::string agent_class_;
            int type_;
            
            Agent *agent_;
            Genome *genome_;
            Chromosome *chromo_;
            std::list< ChromosomeElement* > *chr_;
            AgentTag me_, mother_;
            
            double cp_gene_, rm_gene_;
            double cp_tp_, rm_tp_, dsb_, rm_ltr_;
            
            std::vector< Agent* > population_;
            std::vector< Location > locations_;
    };

inline bool PopulationReader::sawErrors() const
{ return sawErrors_; }

inline bool PopulationReader::done() const
{ return done_; }

inline void PopulationReader::type( int t )
{ type_ = t; }
}
#endif
