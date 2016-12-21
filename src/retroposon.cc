//
// Implementation of a special downstream region, the Retroposon.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "retroposon.hh"

template<> fluke::ObjectCache< fluke::Retroposon >* 
fluke::ObjectCache< fluke::Retroposon >::instance_ = 0;


fluke::ChromosomeElement* 
fluke::Retroposon::clone() const {
    //return new Retroposon( *this );
    Retroposon *rp = ObjectCache< Retroposon >::instance()->borrowObject();
    rp->copy( *this );
    return rp;
}

void 
fluke::Retroposon::copy( const ChromosomeElement &ce ) {
    const Retroposon *tp = dynamic_cast< const Retroposon * >( &ce );
    Downstream::copy( *tp );
}

void
fluke::Retroposon::toPool() {
    ObjectCache< Retroposon >::instance()->returnObject( this );
}

std::string 
fluke::Retroposon::asXmlString() const {
    return "<tposon id=\"" + boost::lexical_cast< std::string >( tag_ ) +
           "\"/>";
}
