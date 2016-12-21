//
// Chromosome of an agent, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_CHROMOSOME_H_
#define _FLUKE_CHROMOSOME_H_

#include "defs.hh"
#include "chromelement.hh"
#include "genome.hh"
#include "bsite.hh"
#include "repeat.hh"
#include "transfac.hh"
#include "retroposon.hh"
#include "ordinary_dstream.hh"
#include "module_dstream.hh"
#include "distribution.hh"
#include "mutate_rates.hh"
// reviewer 2
#include "centromere.hh"

namespace fluke {
   
    /// \class Chromosome
    /// \brief A sequence of b-sites, downstreams, repeats and retroposons.
    ///
    /// The heart of the genome model. Basically, it is a sequence of 
    /// b-sites, downstreams, long terminal repeats and retroposons in any
    /// order. 
    ///
    /// On the chromosome several mutational processes are defined: bsites
    /// appear spontaneously at any location, they are copied around the
    /// chromosome (hypothetical, no mechanism known for in yeast) or get
    /// deleted, retrotransposons copy within the chromosome... FIX
    ///
    /// From this sequence, a transcription network can be built.
    class Chromosome : public CachedElement {
        friend class FixedMutateRates;
        friend class LinearMutateRates;
        friend class UniformMutateRates;
        friend class IterativeMutateRates;
        
        public:
        enum mut_event { DSB, CP_G, RM_G, CP_RP, RM_RP, RM_LTR };
        static const uint NR_RATES = 12;
        
        public:
        /// Chromosome iterator
        typedef std::list< ChromosomeElement* >::iterator ce_iter;
        /// Chromosome const iterator
        typedef std::list< ChromosomeElement* >::const_iterator 
            const_ce_iter;
        /// Chromosome reverse iterator
        typedef std::list< ChromosomeElement* >::reverse_iterator ce_riter;
        /// Chromosome const reverse iterator
        typedef std::list< ChromosomeElement* >::const_reverse_iterator 
            const_ce_riter;
        /// Chromosome tag list/vector/map
        typedef std::vector< uint > tag_container;
        /// Chromosome tag list/vector/map iterator
        typedef tag_container::iterator tag_iter;

        public:
        /// Constructor
        Chromosome();
        /// Constructor, the list is not copied. Memory responsabilities
        /// are moved to \c this.
        Chromosome( Genome *, std::list< ChromosomeElement* > * );
        /// Copy constructor
        explicit Chromosome( const Chromosome & );
        /// Destructor
        ~Chromosome();
        /// Copy a chromosome into \c this
        void copy( const Chromosome & );
        /// Clone a chromosome
        Chromosome* clone() const;
        /// Back to the pool
        virtual void toPool();

        /// All the mutations that can be handled within the chromosome
        /// are performed by invoking this method.
        int mutate();
        /// New, copy and deletion events of binding sites. The expected
        /// number of events is proportional to the abundance of b-sites 
        /// in the chromosome.
        ce_iter bsiteMutate( ce_iter );
        /// Copy and deletion events on the scale of genes. As genes are 
        /// not explicitly in the genome, a gene is defined as an ordinary
        /// downstream or transcription factor element with its sequence of 
        /// directly preceding binding sites.
        ce_iter geneMutate( ce_iter );
        
        /// Retrotransposons are copied around the genome (together with
        /// their accompanying LTRs) via reverse transcriptase. Removal of
        /// retroposons is done via reciprocal recombination.
        ce_iter retroposonMutate( ce_iter );
        /// Long Terminal Repeats are subject to double strand breaks.
        ce_iter repeatMutate( ce_iter );
        /// Retrotransposons may arise newly at rare occasions.
        void newRetrotransposon();
        
        /// Return an iterator to a random element in this chromosome
        boost::tuple< Chromosome*, ce_iter > randChromosomeElement();
        /// Return an iterator to a random element in the genome
        boost::tuple< Chromosome*, ce_iter > randGenomeElement();
        /// Return all the segments
        std::list< Chromosome* > segments();
        /// Append a chromosome to the end of \c this
        void append( Chromosome* );
        
        /// Reset the state of all chromosome elements. It should be 
        /// performed before a new simulation update on a genome
        void reset();
        /// Empty the chromosome, the chromosome elements are not freed
        void clear();
        /// Make sure all flags of caching behaviour are set to update cache
        void recache();
        
        /// Set the parent (genome) of \c this
        void parent( Genome * );
        /// Get the parent of \c this
        Genome* parent() const;

        /// Get a reference to the contents of the chromosome
        const std::list< ChromosomeElement* > & elements() const;

        /// Get the downstream tags present in the chromosome
        tag_container essentialTags() const;
        /// Get the downstream tags of a module
        tag_container moduleTags( int ) const;
        
        /// Get nr retroposons
        int nrRetroposons() const;
        /// Get nr of repeats
        int nrRepeats() const;
        /// Get nr of dsbs
        int nrDoubleStrandBreaks() const;
        /// Get nr of gene cp/rm
        tag_container nrMutations() const;
        /// Get the size of the chromosome (number of elements)
        int size() const;
        /// Is the chromosome empty?
        bool empty() const;
        /// Reviewer 2: has it got one and only one centromere?
        bool oneCentromere() const;
        
        /// Write a (xml) representation to string
        void write( std::ostream & ) const;

        public:
        /// Set copy rate of a retroposon
        void copyRetroposonRate( double );
        /// Get copy rate of a retroposon
        double copyRetroposonRate() const;
        /// Set removal rate of a retroposon
        void removeRetroposonRate( double );
        /// Get removal rate of a retroposon
        double removeRetroposonRate() const;
        /// Set rate of repeat removal
        void removeRepeatRate( double );
        /// Get rate of repeat removal
        double removeRepeatRate() const;
        /// Set DSB recombination rate
        void recombinationRate( double );
        /// Get DSB recombination rate
        double recombinationRate() const;
        /// Set new retroposon rate
        void newRetroposonRate( double );
        /// Get the new retroposon rate
        double newRetroposonRate();
        
        /// Set new binding site generation rate
        void newBsiteRate( double );
        /// Get new binding site generation rate
        double newBsiteRate() const;
        /// Set binding site deletion rate
        void removeBsiteRate( double );
        /// Get binding site deletion rate
        double removeBsiteRate() const;
        /// Set binding site copy rate
        void copyBsiteRate( double );
        /// Get binding site copy rate
        double copyBsiteRate() const;
        
        /// Set gene (up- and downstream) copy rate
        void copyGeneRate( double );
        /// Get gene copy rate
        double copyGeneRate() const;
        /// Set gene (up- and downstream) deletion rate
        void removeGeneRate( double );
        /// Get gene deletion rate
        double removeGeneRate() const;
        
        public:
        /// Set mutational step
        void mutationStep( double );
        /// Get mutational step
        double mutationStep() const;
        /// Set DSB recombination step
        void dsbStep( double );
        /// Get DSB recombination step
        double dsbStep() const;
        /// Set retroposon step
        void retroStep( double );
        /// Get retroposon step
        double retroStep() const;
        /// Set mutation rate of rates
        void mutationRate( double );
        /// Get mutation rate of rates
        double mutationRate() const;

        public:
        /// Set mutational scheme
        static void mutationScheme( MutateRates* );
        
        public:
        class IsRetroposon : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsRetroposon() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Retroposon ); }
        };
        
        class IsTrueDstream : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsTrueDstream() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( TranscriptionFactor ) or
                     typeid( *ce ) == typeid( OrdinaryDownstream ) or
                     typeid( *ce ) == typeid( ModuleDownstream ); }
        };
        
        class IsBindingSite : 
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsBindingSite() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( BindingSite ); }
        };
        
        class IsRepeat :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsRepeat() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Repeat ); }
        };
        
        class IsDoubleStrandBreak :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsDoubleStrandBreak() {}
            bool operator()( ChromosomeElement *ce ) const {
                Repeat *aux = dynamic_cast< Repeat* >( ce );
                if( aux ) {
                    return aux->hasDSB();
                } else {
                    return false;
                }
            }
        };

        class IsCentromere :
            public std::unary_function< ChromosomeElement*, bool > {
            public:
            IsCentromere() {}
            bool operator()( ChromosomeElement *ce ) const
            { return typeid( *ce ) == typeid( Centromere ); }
        };
        
        private:
        ce_iter upstreamSelect( ce_iter );
        int nrRetroposons( ce_iter, ce_iter ) const;
        void copyRates( const Chromosome &, Chromosome & ) const;
        // overloading list methods coz of length caching; insert, splice
        ce_iter insert( ce_iter, ChromosomeElement* );
        void splice( ce_iter, std::list< ChromosomeElement* >, uint, uint );
        
        private:
        Genome *parent_;
        std::list< ChromosomeElement* > *chro_;
        std::vector< uint > mut_events_;
        mutable uint nr_retroposons_, nr_ltr_, len_;
        mutable bool update_retro_, update_ltr_, update_len_;
            
        // mutation rates 
        double cp_tp_rate_;
        double rm_tp_rate_;
        double rm_ltr_rate_;
        double new_tp_rate_;

        double nw_bs_rate_;
        double cp_bs_rate_;
        double rm_bs_rate_;

        double cp_gene_rate_;
        double rm_gene_rate_;
        
        // rate of a dsb occurrence * rate of recombination repair
        double dsb_recombination_;
        
        // rates to evolve mutation rates
        double mut_step_;
        double dsb_step_;
        double retro_step_;
        double mut_rate_;
        
        private:
        static MutateRates *rate_mutator_;
    };

    /// Overloaded \c << operator for easy writing to streams.
    inline std::ostream& operator<<( std::ostream& os, const Chromosome& chro )
    { chro.write( os ); return os; }

    inline void Chromosome::parent( Genome *p )
    { parent_ = p; }

    inline Genome* Chromosome::parent() const 
    { return parent_; }

    inline boost::tuple< fluke::Chromosome*, fluke::Chromosome::ce_iter >
    Chromosome::randGenomeElement()
    { return parent_->randElement(); }

    inline const std::list< ChromosomeElement* >& Chromosome::elements() const
    { return *chro_; }

    inline int Chromosome::nrDoubleStrandBreaks() const 
    { return mut_events_[ DSB ]; }

    inline std::vector< uint > Chromosome::nrMutations() const
    { return mut_events_; }
    
    inline bool Chromosome::empty() const 
    { return chro_->empty(); }
    
    inline void Chromosome::clear()
    { chro_->clear(); }
}
#endif

