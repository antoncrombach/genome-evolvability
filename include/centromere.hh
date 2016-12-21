//
// Reviewer 2 wants a centromere...
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_CENTROMERE_H_
#define _FLUKE_CENTROMERE_H_

#include "defs.hh"
#include "pool.hh"
#include "chromelement.hh"

namespace fluke {

    /// \class Centromere
    /// \brief Centromere of a chromosome
    ///
    /// Centromeres, these poorly understood bastards. 
    class Centromere : public ChromosomeElement {
        public:
        /// Constructor.
        Centromere();
        /// Copy constructor.
        explicit Centromere( const Centromere & );
        /// Destructor.
        virtual ~Centromere() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another Centromere into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return Centromere to pool
        virtual void toPool();
        
        /// Dummy mutate
        virtual int mutate();
       
        /// Write a html representation.
        virtual std::string asXmlString() const;
        /// Write a string representation
        virtual std::string asString() const;
    };

    inline std::string Centromere::asString() const
    { return "centro"; }

    inline int Centromere::mutate() 
    { return 0; }
    
}
#endif
