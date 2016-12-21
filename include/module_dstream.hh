//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_MODULE_DSTREAM_H_
#define _FLUKE_MODULE_DSTREAM_H_

#include "defs.hh"
#include "downstream.hh"
#include "shortseq.hh"

namespace fluke {

    /// \class ModuleDownstream
    /// \brief Module downstreams are genes that can be categorised in
    /// different modules.
    class ModuleDownstream : public Downstream {
        public:
        /// Constructor.
        ModuleDownstream() : Downstream(), module_( -1 ) {}
        /// Constructor with module tag and identification tag
        ModuleDownstream( int t, int m ) : Downstream( t ), module_( m ) {}
        /// Copy constructor.
        explicit ModuleDownstream( const ModuleDownstream &md )
        { ModuleDownstream::copy( md ); }
        /// Destructor.
        virtual ~ModuleDownstream() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another ModuleDownstream into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return to pool.
        virtual void toPool();
        
        /// Mutate the protein(?).
        virtual int mutate();
        /// Set the module of the downstream.
        void module( int );
        /// Get the module of the downstream.
        int module() const;

        /// Write a html representation.
        virtual std::string asXmlString() const;
        /// Write a string representation
        virtual std::string asString() const;
        /// Write a label for graphs.
        virtual std::string writeLabel() const;

        private:
        int module_;        
    };
    
    inline std::string ModuleDownstream::writeLabel() const
    { return boost::lexical_cast< std::string>( module_ ); }
    
    inline std::string ModuleDownstream::asString() const
    { return boost::lexical_cast< std::string >( tag_ ) + ":" +
        boost::lexical_cast< std::string>( module_ ); }

    inline void ModuleDownstream::module( int m )
    { module_ = m; }
    
    inline int ModuleDownstream::module() const
    { return module_; }
}
#endif

