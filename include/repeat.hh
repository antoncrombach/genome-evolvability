//
// Chromosome building block, part of the genome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_REPEAT_H_
#define _FLUKE_REPEAT_H_

#include "defs.hh"
#include "pool.hh"
#include "chromelement.hh"

namespace fluke {

    /// \class Repeat
    /// \brief Long Terminal Repeat of a retrotransposon.
    ///
    /// Repeats are a known mechanism for the occurrence of double-strand
    /// breaks (DSBs). 
    class Repeat : public ChromosomeElement {
        public:
        /// Constructor.
        Repeat();
        /// Copy constructor.
        explicit Repeat( const Repeat & );
        /// Destructor.
        virtual ~Repeat() {}
        /// Clone \c this.
        virtual ChromosomeElement* clone() const;
        /// Copy another Repeat into \c this.
        virtual void copy( const ChromosomeElement & );
        /// Return Repeat to pool
        virtual void toPool();
        
        /// Dummy mutate or incorporate possibility of DSB?
        virtual int mutate();

        /// Fix any double-strand break.
        void repairDSB();
        /// Induce a double-strand break.
        void induceDSB();
        /// Returns if it has a double-strand break.
        bool hasDSB() const;
        
        /// Write a html representation.
        virtual std::string asXmlString() const;
        /// Write a string representation
        virtual std::string asString() const;
        
        private:
        bool dsb_;
    };

    inline void Repeat::induceDSB()
    { dsb_ = true; }    

    inline void Repeat::repairDSB() 
    { dsb_ = false; }

    inline bool Repeat::hasDSB() const 
    { return dsb_; }

    inline std::string Repeat::asString() const
    { return "ltr"; }

    inline int Repeat::mutate() 
    { return 0; }
}
#endif

