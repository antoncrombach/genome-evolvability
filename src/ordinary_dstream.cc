//
// Implementation of an ordinary downstream region.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "ordinary_dstream.hh"


template<> fluke::ObjectCache< fluke::OrdinaryDownstream >* 
fluke::ObjectCache< fluke::OrdinaryDownstream >::instance_ = 0;


fluke::ChromosomeElement* 
fluke::OrdinaryDownstream::clone() const {
    //return new OrdinaryDownstream( *this );
    OrdinaryDownstream *tp =
        ObjectCache< OrdinaryDownstream >::instance()->borrowObject();
    tp->copy( *this );
    return tp;
}

void
fluke::OrdinaryDownstream::toPool() {
    ObjectCache< OrdinaryDownstream >::instance()->returnObject( this );
}

std::string 
fluke::OrdinaryDownstream::asXmlString() const {
    return "<dstream id=\"" + boost::lexical_cast< std::string >( tag_ ) + 
           "\"/>\n";
}
