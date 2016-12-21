//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_RETROPOSON_H_
#define _FLUKE_RETROPOSON_H_

#include "defs.hh"
#include "downstream.hh"
#include "shortseq.hh"

namespace fluke {

    /// \class Retroposon
    /// \brief Retrotransposable element.
    ///
    /// Retroposons are a known mechanism for gene duplication and for the 
    /// occurrence of double-strand breaks (DSBs). 
    /// A retroRetroposon is always flanked by long terminal repeats, 
    /// represented by elements of \c Repeat.
    class Retroposon : public Downstream {
        public:
        /// Constructor.
        Retroposon() : Downstream() {}
        /// Constructor. The \c int is a unique identifier.
        Retroposon( uint t ) : Downstream( t ) {}
        /// Copy constructor.
        Retroposon( const Retroposon &rp )
        { Retroposon::copy( rp ); }
        /// Destructor.
        virtual ~Retroposon() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another retroposon into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual void toPool();
        
        /// Mutate the protein.
        virtual int mutate();

        /// Write a html representation.
        virtual std::string asXmlString() const;
        /// Write a string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string writeLabel() const;
    };

    inline std::string Retroposon::writeLabel() const
    { return std::string( "xx" ); }
    
    inline std::string Retroposon::asString() const
    { return boost::lexical_cast< std::string >( tag_ ); }

    inline int Retroposon::mutate() 
    { return 0; }
}
#endif

