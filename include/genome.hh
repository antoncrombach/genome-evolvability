//
// Genome of an agent.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_GENOME_H_
#define _FLUKE_GENOME_H_

#include "defs.hh"

namespace fluke {

    /// \class Genome
    /// \brief Container of chromosomes.
    ///
    /// The genome is one of the three levels in a cell that we distinguish. It 
    /// harbours a high level object model of a genome. The genome consists of
    /// chromosome (which consist of genes etc) and has several mutational 
    /// operators defined on them. On chromosome level the mutational process
    /// is chromosomal tail swapping. Such mutations may occur if double
    /// stranded breaks are not repaired correctly.
    class Genome {
        public:
            /// Chromosome iterator 
            typedef std::list< Chromosome* >::iterator chromos_iter;
            /// Reverse chromosome iterator
            typedef std::list< Chromosome* >::reverse_iterator chromos_riter;
            /// Constant chromosome iterator
            typedef std::list< Chromosome* >::const_iterator 
                const_chromos_iter;
            /// Constant reverse chromosome iterator
            typedef std::list< Chromosome* >::const_reverse_iterator 
                const_chromos_riter;

        public:
            /// (Dummy) constructor
            Genome();
            /// Constructor with list of chromosomes
            Genome( std::list< Chromosome* > * );
            /// Copy constructor
            Genome( const Genome & );
            /// Get clone of this genome (OO pattern \a prototype)
            Genome* clone() const;
            /// Copy given genome into \c this
            void copy( const Genome & );
            /// Destructor
            ~Genome();

            /// Duplicate the genome
            void duplicate();
            /// Split the genome in two
            Genome* split();
            /// Mutate the genome. Returns number of mutations that occurred
            int mutate(); 
            /// Write the genome to an output stream (in xml format)
            void write( std::ostream & ) const;

            /// Return all the tags present in the genome (from downstreams)
            std::vector< uint > essentialTags() const;
            /// Return all the tags present of a certain downstream module
            std::vector< uint > moduleTags( int ) const;

            /// Return a pointer to any element but a repeat or retroposon
            boost::tuple< fluke::Chromosome*, 
                std::list< ChromosomeElement* >::iterator > randElement();
            
            /// Get number of chromosomes
            int size() const;
            /// Get total number of chromosome elements (i.e. binding sites,
            /// transposons, module downstream, ordinary downstreams)
            int fullSize() const;
            /// Is the genome empty?
            bool empty() const;
            /// Are there any mutations worth mentioning
            bool hasMutation() const;

            /// Get constant reference to the chromosomes. Used by observers.
            const std::list< Chromosome* > & chromosomes() const;
            /// Get nr of transposons
            int nrRetroposons() const;
            /// Get nr of single ltrs
            int nrRepeats() const;
            /// Get nr of double strand breaks
            int nrDoubleStrandBreaks() const;
            /// Get nr of dsbs parent had in diploid phase
            boost::tuple< int, int > nrDsbParent() const;
            /// Get nr of dsbs, gene cp/rm
            std::vector< uint > nrMutations() const;
            /// Reviewer 2: one centromere is healthy
            bool oneCentromere() const;

        private:
            void clear();
            
        private:
            std::list< Chromosome* > *chromos_;
            // The idea of these two vectors is that they are not changed during
            // the life of an individual, and only get new values at mutating
            // (=reproduction)
            std::vector< uint > nr_dsbs_parent_;
            std::vector< uint > nr_mutations_;
    };
    
    /// Overloaded \c << operator for easy writing to streams.
    inline std::ostream& operator<<( std::ostream& os, const Genome& genome )
    { genome.write( os ); return os; }

    inline int Genome::size() const
    { return chromos_->size(); }

    inline bool Genome::empty() const
    { return chromos_->empty(); }

    inline const std::list< fluke::Chromosome* > & Genome::chromosomes() const 
    { return *chromos_; }
    
    inline boost::tuple< int, int > Genome::nrDsbParent() const
    { return boost::make_tuple( nr_dsbs_parent_[ 0 ], nr_dsbs_parent_[ 1 ] ); }
    
    inline std::vector< uint > Genome::nrMutations() const
    { return nr_mutations_; }
}
#endif

