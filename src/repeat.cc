//
// Implementation of a special chromosome element, the repeat.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "repeat.hh"
#include "pool.hh"

template<> fluke::ObjectCache< fluke::Repeat >* 
fluke::ObjectCache< fluke::Repeat >::instance_ = 0;

fluke::Repeat::Repeat() : ChromosomeElement(), dsb_( false ) {}

fluke::Repeat::Repeat( const Repeat &tp ) : ChromosomeElement( tp ) {
    Repeat::copy( tp );
}

fluke::ChromosomeElement* 
fluke::Repeat::clone() const {
    //return new Repeat( *this );
    Repeat *ltr = ObjectCache< Repeat >::instance()->borrowObject();
    ltr->dsb_ = dsb_;
    return ltr;
}

void 
fluke::Repeat::copy( const ChromosomeElement &ce ) {
    const Repeat *ltr = dynamic_cast< const Repeat * >( &ce );
    dsb_ = ltr->dsb_;
}

void
fluke::Repeat::toPool() {
    ObjectCache< Repeat >::instance()->returnObject( this );
}

std::string 
fluke::Repeat::asXmlString() const { 
    std::string aux( "<repeat dsb=" );
    aux += ( dsb_? "\"yes\"": "\"no\"" );
    aux += "/>\n";
    return aux;
}

