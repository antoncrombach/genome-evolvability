//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_ORDINARY_DSTREAM_H_
#define _FLUKE_ORDINARY_DSTREAM_H_

#include "defs.hh"
#include "downstream.hh"
#include "shortseq.hh"

namespace fluke {

    /// \class OrdinaryDownstream
    /// \brief Ordindary downstreams are all the genes not involved in 
    /// regulation. 
    class OrdinaryDownstream : public Downstream {
        public:
        /// Constructor.
        OrdinaryDownstream() : Downstream() {}
        /// Constructor with tag
        OrdinaryDownstream( int t ) : Downstream( t ) {}
        /// Copy constructor.
        OrdinaryDownstream( const OrdinaryDownstream &tp ) 
        { OrdinaryDownstream::copy( tp ); }
        /// Destructor.
        virtual ~OrdinaryDownstream() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another OrdinaryDownstream into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual void toPool();
        
        /// Mutate the protein(?).
        virtual int mutate();

        /// Write a html representation.
        virtual std::string asXmlString() const;
        /// Write string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string writeLabel() const;
    };
    
    inline void OrdinaryDownstream::copy( const ChromosomeElement &ce ) 
    { Downstream::copy( ce ); }

    inline std::string OrdinaryDownstream::asString() const
    { return boost::lexical_cast< std::string >( tag_ ); }
    
    inline std::string OrdinaryDownstream::writeLabel() const
    { return std::string( "G" ); }

    inline int OrdinaryDownstream::mutate()
    { return 0; }
}
#endif

