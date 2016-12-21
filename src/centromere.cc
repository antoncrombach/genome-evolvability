//
// Implementation of the centromere
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "centromere.hh"
#include "pool.hh"

template<> fluke::ObjectCache< fluke::Centromere >* 
fluke::ObjectCache< fluke::Centromere >::instance_ = 0;

fluke::Centromere::Centromere() : ChromosomeElement() {}

fluke::Centromere::Centromere( const Centromere &tp ) 
    : ChromosomeElement( tp ) {
    Centromere::copy( tp );
}

fluke::ChromosomeElement* 
fluke::Centromere::clone() const {
    //return new Centromere( *this );
    Centromere *ltr = ObjectCache< Centromere >::instance()->borrowObject();
    return ltr;
}

void 
fluke::Centromere::copy( const ChromosomeElement &ce ) {}

void
fluke::Centromere::toPool() {
    ObjectCache< Centromere >::instance()->returnObject( this );
}

std::string 
fluke::Centromere::asXmlString() const { 
    return std::string( "<centromere/>\n" );
}

