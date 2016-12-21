//
// Interface of a genome agent with gene modules.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_MODULEAGENT_H_
#define _FLUKE_MODULEAGENT_H_

#include "defs.hh"
#include "agent.hh"
#include "genome.hh"
#include "chromosome.hh"


namespace fluke {

    /// \class ModuleAgent
    /// \brief Agent with a genome and gene modules.
    ///
    /// Our first \a real simulations are performed with this agent. It has
    /// a genome consisting of genes. The genes can be subdivided in essential
    /// genes and modules. The modules are genes with a certain environmental
    /// benefit.
    class ModuleAgent : public Agent {
        friend class AgentReader;
        public:
        typedef std::vector< Chromosome::tag_container >::iterator module_iter;
        typedef std::vector< Chromosome::tag_container >::const_iterator
            const_module_iter;
        
        public:
        /// Constructor
        ModuleAgent( int, Genome* );
        /// Copy constructor
        ModuleAgent( const ModuleAgent & );
        /// Destructor
        virtual ~ModuleAgent();
        /// Cloning the ModuleAgent
        virtual Agent* clone() const;
        /// Copying another ModuleAgent into \c this
        virtual void copy( const Agent & );
        
        /// Dummy
        virtual void initialise();
        /// Perform an update step (one simulation timestep)
        virtual void step( Population & );
        /// Reproduction
        virtual Agent* sibling();
        
        /// Get the score of the agent
        virtual double score() const;
        /// Calculate the score of the individual
        void evaluate( const Environment & );
        /// Count all genes
        void countGenes();
        /// Count essential genes
        void countEssentialGenes();
        /// Count genes of all modules
        void countModuleGenes();
        
        /// Write a (xml) representation of the agent to stream
        virtual void write( std::ostream & ) const;
        
        /// Get the genome
        const Genome& genome() const;
        /// Get more gene counts
        const Chromosome::tag_container nrEssentialGenes();
        /// Get gene counts
        const std::vector< Chromosome::tag_container > nrModuleGenes();
        
        /// Get nr modules
        int nrModules() const;
        /// Get nr transposons
        int nrRetroposons() const;
        /// Get nr repeats
        int nrRepeats() const;
        /// Get nr of dsbs
        int nrDoubleStrandBreaks() const;
        /// Get nr of dsbs in diploid parent
        boost::tuple< int, int > nrDsbParent() const;
        /// Get a few mutation counters
        std::vector< uint > nrMutations() const;
        /// Get genotypical distance to target
        int distance() const;
        /// Get distance of parent
        int distanceParent() const;
        /// Get size of agent's genome
        int size() const;
        /// Get size of parent's genome
        int sizeParent() const;
        
        public:
        /// Set birth rate of agent
        static void birthRate( float );
        /// Get birth rate of agent
        static float birthRate();
        /// Set death rate of agent
        static void deathRate( float );
        /// Get death rate of agent
        static float deathRate();
        /// Fill a tag reference vector
        static void referenceTags( const ModuleAgent & );
        /// Fill a lookup table
        static void lookupTable();
        /// Set maximum distance
        static void maxDistance( int );
        /// Get maximum distance
        static int maxDistance();
        /// Set max genome size, after which a penalty is applied
        static void maxGenomeSize( int );
        /// Set max # retroposons, after which a penalty is applied
        static void maxTposons( int );
        /// Set penalty per unit over the max size/number
        static void genomePenaltyRate( double );
        static void retroposonPenaltyRate( double );
        
        protected:
        /// Calculate the score of the essential genes
        int essentialsScore( const Environment &env );
        /// Calculate the score of the module genes
        int modulesScore( const Environment &env );
        /// Give it a penalty if too big
        double penalty( double ) const;
            
        protected:
        /// Current score (a.k.a. fitness)
        int distance_;
        /// Score of parent
        int distance_parent_;
        /// Size of parent genome
        int size_parent_;
        /// Genome of the agent
        Genome *genome_;
        /// Flag
        bool inventorised_;
        /// Container with nr of module gene tags of \c this genome
        std::vector< Chromosome::tag_container > mod_tags_now_;
        /// Container with nr of essential gene tags of \c this genome
        Chromosome::tag_container ess_tags_now_;
            
        protected:
        /// Birth rate
        static float birth_rate_;
        /// Death rate
        static float death_rate_;
        /// Reference container with module gene tags
        static std::vector< Chromosome::tag_container > module_tags_;
        /// Reference container with essential gene tags
        static Chromosome::tag_container essential_tags_;
        /// Maximum distance we look at
        static int max_dist_;
        /// Coefficient for max distance
        static double dist_coeff_;
        /// Penalty if genome size is greater
        static int penalty_genome_;
        /// Penalty if more retroposons than allowed
        static int penalty_tposons_;
        /// Penalty coefficient for genome
        static double penalty_genome_rate_;
        /// Penalty coefficient for retroposon
        static double penalty_tposons_rate_;
    };
    
    inline const Genome & ModuleAgent::genome() const
    { return *genome_; }
    
    inline int ModuleAgent::nrModules() const
    { return module_tags_.size(); }

    inline const Chromosome::tag_container ModuleAgent::nrEssentialGenes()
    { if( !inventorised_ ) countGenes(); return ess_tags_now_; }

    inline const std::vector< Chromosome::tag_container >
    ModuleAgent::nrModuleGenes()
    { if( !inventorised_ ) countGenes(); return mod_tags_now_; }
    
    inline int ModuleAgent::nrRetroposons() const
    { return genome_->nrRetroposons(); }

    inline int ModuleAgent::nrRepeats() const
    { return genome_->nrRepeats(); }

    inline int ModuleAgent::nrDoubleStrandBreaks() const
    { return genome_->nrDoubleStrandBreaks(); }
    
    inline boost::tuple< int, int > ModuleAgent::nrDsbParent() const
    { return genome_->nrDsbParent(); }
    
    inline std::vector< uint > ModuleAgent::nrMutations() const
    { return genome_->nrMutations(); }
    
    inline int ModuleAgent::distance() const
    { return distance_; }

    inline int ModuleAgent::distanceParent() const
    { return distance_parent_; }

    inline int ModuleAgent::size() const
    { return genome_->fullSize(); }
    
    inline int ModuleAgent::sizeParent() const
    { return size_parent_; }
}
#endif

