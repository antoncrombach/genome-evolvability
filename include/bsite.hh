//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_BSITE_H_
#define _FLUKE_BSITE_H_

#include "defs.hh"
#include "chromelement.hh"
#include "shortseq.hh"

namespace fluke {

    /// \class BindingSite
    /// \brief Chromosome element, represents a transcription factor binding 
    /// site (TFBS).
    ///
    /// A binding site (b-site) is a short sequence in the upstream region of 
    /// a gene, about 5 to 12 nucleotides. Transcription factors may bind to a
    /// site and activate or repress the transcription of the downstream 
    /// sequence (the gene).
    class BindingSite : public ChromosomeElement {
        public:
        /// Constructor
        BindingSite();
        /// Constructor, give the short sequence in integer encoding
        BindingSite( label );
        /// Copy constructor
        explicit BindingSite( const BindingSite & );
        /// Destructor
        virtual ~BindingSite();
        /// Clone the b-site
        virtual ChromosomeElement* clone() const;
        /// Copy another b-site into \c this
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual void toPool();
        
        /// Point mutate the short sequence
        virtual int mutate();
        /// Return integer encoding the nucleotide short sequence
        label tfbs() const;
        /// Set the nucleotide short sequence
        void tfbs( label );

        /// Write HTML representation of b-site
        virtual std::string asXmlString() const;
        /// Write string representation
        virtual std::string asString() const;

        public:
        /// Set short sequence manager
        static void shortSeqManager( ShortSeqManager * );
        /// Get short sequence manager
        static ShortSeqManager* shortSeqManager();
            
        private:
        label tfbs_;

        private:
        static ShortSeqManager *ssm_;
    };

    inline label BindingSite::tfbs() const
    { return tfbs_; }

    inline void BindingSite::tfbs( label l )
    { tfbs_ = l; }

}
#endif

