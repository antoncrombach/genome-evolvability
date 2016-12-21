//
// Factory for creating a genome and other model entities.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_FACTORY_H_
#define _FLUKE_FACTORY_H_

#include "defs.hh"
#include "config.hh"

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/XMLString.hpp>


namespace fluke {

    XERCES_CPP_NAMESPACE_USE
    
    /// \class Factory
    /// \brief Model parts are built in the \c Factory.
    ///
    /// Given the configuration, parts of the model can be built in the factory.
    class Factory {
        friend class AgentReader;
        friend class PopulationReader;
        public:
            /// (Dummy) constructor
            Factory();
            /// Constructor with configuration
            Factory( Config *cg );
            /// Destructor
            ~Factory();

            /// Create binding site. A binding site manager is created if 
            /// necessary.
            BindingSite* bsite();
            /// Create transcription factor. A transcription factor manager is
            /// created if necessary.
            TranscriptionFactor* transfac();
            /// Create downstream with module identification.
            ModuleDownstream* moduleDstream( int );
            /// Create an ordinary downstream.
            OrdinaryDownstream* ordinaryDstream();
            /// Create a transposon (a specialised type of downstream).
            Retroposon* retroposon();
            /// Create a repeat (long terminal repeat)
            Repeat* repeat();
            /// Reviewer 2 wanted a centromere
            Centromere* centromere();
            
            /// Create a chromosome. It is filled with binding sites,
            /// downstreams and transposons (not interrupting any upstream)
            /// according to the settings in the configuration
            Chromosome* chromosome( int );
            /// Create a genome
            Genome* genome( int );
            /// Create mutational scheme we want to use
            MutateRates* mutateRates( int);            

            /// Create an agent. The type of agent is given in the 
            /// configuration
            Agent* agent( int );
            /// Create the population, still needs to be initialised.
            Population* population();
            /// Create the environment, still needs to be initialised.
            /// The type is given in the configuration.
            Environment* environment();

            /// Create a binding site manager once, return that one instance
            /// every time the method is called (singleton behaviour)
            ShortSeqManager* bsiteManager();
            /// Create a transcription factor manager once, return that one
            /// instance every time the method is called (singleton behaviour)
            ShortSeqManager* transfacManager();

            /// Create a scaling scheme used to scale raw scores of the agents
            ScalingScheme* scalingScheme();
            /// Create a selection scheme used to select an agent for
            /// reproduction based on its score.
            SelectionScheme* selectionScheme();

            /// Create an observer manager. It keeps track of observers and to
            /// which subject they are linked.
            ObserverManager* observerManager();

        private:
            Agent* readAgent( std::string, int );
            boost::tuple< std::vector< Agent* >, std::vector< Location > > 
                readPopulation( std::string, int );
            void readAgentConfigurations();
            void setConfiguration( int );
            
            // build correct # of elements according to config
            std::vector< std::list< ChromosomeElement* > > *
                chromosomeParts( int );
            // organise chromosome up to certain degree
            std::list< ChromosomeElement* >* organiseChromosome( 
                std::vector< std::list< ChromosomeElement* > > *, double );
            
            void resetTag();
            int tag();

        private:
            ShortSeqManager *bs_mng_, *tf_mng_;
            Config *conf_;
            int downstream_tag_;
    };

    inline void Factory::resetTag()
    { downstream_tag_ = 0; }
    
    inline int Factory::tag()
    { return downstream_tag_++; }
}
#endif

