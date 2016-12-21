//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_TRANSCRIPTION_FACTOR_H_
#define _FLUKE_TRANSCRIPTION_FACTOR_H_

#include "defs.hh"
#include "downstream.hh"
#include "shortseq.hh"
#include "bsite.hh"

namespace fluke {

    /// \class TranscriptionFactor
    /// \brief Transcription factors are all the genes that may be involved in 
    /// regulation. 
    class TranscriptionFactor : public Downstream {
        public:
        /// Constructor.
        TranscriptionFactor() : Downstream(), tf_( -1 ) {}
        /// Constructor with tag and label to give.
        TranscriptionFactor( uint t, label l ) : Downstream( t ), tf_( l ) {}
        /// Copy constructor.
        explicit TranscriptionFactor( const TranscriptionFactor & );
        /// Destructor.
        virtual ~TranscriptionFactor();
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another TranscriptionFactor into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual void toPool();
        
        /// Mutate the protein; its recognition site.
        virtual int mutate();
        /// Does \c this bind to the given binding site.
        bool binds( const BindingSite & ) const;
        /// Get terminal repeat in integer encoding.
        label transFac() const;
        /// Set terminal repeat in integer encoding.
        void transFac( label );
        
        /// Write a xml representation.
        virtual std::string asXmlString() const;
        /// Write a string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string writeLabel() const;

        public:
        /// Set terminal repeat manager.
        static void shortSeqManager( ShortSeqManager * );

        private:
        label tf_;

        private:
        static ShortSeqManager *ssm_;
    };
    
    inline bool TranscriptionFactor::binds( const BindingSite &bs ) const 
    { return ssm_->similar( tf_, bs.tfbs() ); }
    
    inline label TranscriptionFactor::transFac() const
    { return tf_; }

    inline void TranscriptionFactor::transFac( label l )
    { tf_ = l; }

    inline std::string TranscriptionFactor::asString() const
    { return ssm_->strShortSeq( tf_ ); }
    
    inline std::string TranscriptionFactor::writeLabel() const
    { return ssm_->strShortSeq( tf_ ); }
}
#endif

