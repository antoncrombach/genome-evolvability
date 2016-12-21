//
// Abstract class of downstream regions
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_DOWNSTREAM_H_
#define _FLUKE_DOWNSTREAM_H_

#include "chromelement.hh"

namespace fluke {

    /// \class Downstream
    /// \brief Abstract class of downstream regions
    class Downstream : public ChromosomeElement {
        public:
            /// Virtual destructor.
            virtual ~Downstream() {};
            /// Signature of clone method.
            virtual ChromosomeElement* clone() const = 0;
            /// Signature of copy method.
            virtual void copy( const ChromosomeElement & );
            /// Signature of return to pool method.
            virtual void toPool() = 0;
            
            /// Set tag
            void tag( uint );
            /// Get tag
            uint tag() const;
            
            /// Signature of method for writing a label on a graph vertex.
            virtual std::string writeLabel() const = 0;
            
        protected:
            /// Hidden constructor.
            Downstream() : ChromosomeElement(), tag_( 0 ) {};
            /// Hidden constructor with tag
            Downstream( uint t ) : ChromosomeElement(), tag_( t ) {};
            /// Hidden copy constructor.
            Downstream( const Downstream &ds ) 
                : ChromosomeElement( ds ), tag_( ds.tag_ ) {};
        
        protected:
            /// Identification (e.g. for ancestor tracking)
            uint tag_;
    };

    inline void Downstream::tag( uint t )
    { tag_ = t; }

    inline uint Downstream::tag() const
    { return tag_; }
    
    inline void Downstream::copy( const ChromosomeElement &ce ) 
    { tag_ = dynamic_cast< const Downstream* >( &ce )->tag_; }
}
#endif
