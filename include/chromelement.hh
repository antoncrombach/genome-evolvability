//
// Abstract building block of Chromosome.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_CHROMELEMENT_H_
#define _FLUKE_CHROMELEMENT_H_

#include "defs.hh"

namespace fluke {

    /// \class ChromosomeElement
    /// \brief Abstract element of a chromosome.
    ///
    /// In \c fluke chromosomes consist of \em atoms, called chromosome 
    /// elements. 
    class ChromosomeElement : public CachedElement {
        public:
        /// Virtual destructor.
        virtual ~ChromosomeElement() {};
        /// Signature of clone method.
        virtual ChromosomeElement* clone() const = 0;
        /// Signature of copy method.
        virtual void copy( const ChromosomeElement & ) = 0;
        
        /// Signature of mutate method.
        virtual int mutate() = 0;
        
        /// Deactivate the element.
        void inactivate();
        /// Activate the element.
        void activate();
        /// Is the element \em activated?
        bool isActive() const;
        
        /// Write to an output stream.
        void write( std::ostream & ) const;
        /// Write an XML string.
        virtual std::string asXmlString() const = 0;
        /// write a string representation
        virtual std::string asString() const = 0;

        protected:
        /// Hidden constructor.
        ChromosomeElement() : CachedElement(), active_( true ) {};
        /// Hidden copy constructor.
        explicit ChromosomeElement( const ChromosomeElement &ce ) 
        : CachedElement(), active_( ce.active_ ) {};

        private:
        /// Flag signalling if the element is active.
        bool active_;
    };

    /// Overloaded \c << operator for easy writing to streams.
    inline std::ostream& operator<<( std::ostream& os, 
            const ChromosomeElement& chro )
    { chro.write( os ); return os; }

    inline void ChromosomeElement::inactivate()
    { active_ = false; }

    inline void ChromosomeElement::activate()
    { active_ = true; }

    inline bool ChromosomeElement::isActive() const
    { return active_; }

    inline void ChromosomeElement::write( std::ostream &os ) const 
    { os << asXmlString(); }
}
#endif

