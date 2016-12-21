//
// Short sequence repository, a kind of flyweight manager.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_SHORTSEQ_H_
#define _FLUKE_SHORTSEQ_H_

#include "defs.hh"

namespace fluke {

    //
    // Reference counting policy:
    // - methods expect allocated labels.
    // - methods free the labels themselves.
    //
    /// \class ShortSeqManager
    /// \brief Manages small pieces of DNA
    ///
    /// The applied technique is based on handing out integers as references,
    /// and keeping a repository with all possible pieces of DNA. The references
    /// are indices to the repository. To know how many of each shortseqs are
    /// used, reference counting is performed.
    ///
    /// In the future we may change to using short sequences of differents
    /// lengths, which implies rethinking the inner workings of this class.
    class ShortSeqManager {
        public:
            /// Constructor requiring an alphabet, a length and a point
            /// mutation rate
            ShortSeqManager( std::string, int, double );
            /// Destructor
            ~ShortSeqManager();
        
            /// Allocate the shortseq referenced to by the label 
            void allocShortSeq( label );
            /// Free the shortseq referenced to by the label. Assuming the label
            /// has been allocated in the past
            void freeShortSeq( label );

            /// Allocate random short sequence and return its reference
            label allocRandShortSeq();
            /// Mutate a short seq, return the new sequence and the number of 
            /// point mutations
            boost::tuple< label, int > mutateShortSeq( label );
            /// Are these two short sequences very similar?
            bool similar( label, label );

            /// Get a reference to the reference counts of all the short seqs
            const std::vector< int >& referenceCounts() const;

            /// Return the short sequence, given a reference
            std::string strShortSeq( label ) const;

        public:
            /// Set the maximum hamming distance, such that two short sequences
            /// are still considered similar
            static void maxDistance( int );
            /// Get the maximum hamming distance
            static int maxDistance();

        private:
            void generateShortSeqs();
            void init( std::string, int, int, double );
            
        private:
            std::string alphabet_;
            int length_;
            double mutation_rate_;
            std::vector< std::string > *short_seqs_;
            std::vector< int > *references_;

            randrange_gen_type rand_length_;
            randrange_gen_type rand_repository_;
            randrange_gen_type rand_alphabet_;

        private:
            static int max_hamming_;
    };

}
#endif

