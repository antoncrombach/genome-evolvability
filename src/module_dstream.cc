//
// Implementation of a module downstream region.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "module_dstream.hh"


template<> fluke::ObjectCache< fluke::ModuleDownstream >* 
fluke::ObjectCache< fluke::ModuleDownstream >::instance_ = 0;


fluke::ChromosomeElement* 
fluke::ModuleDownstream::clone() const {
    //return new ModuleDownstream( *this );
    ModuleDownstream *md = 
        ObjectCache< ModuleDownstream >::instance()->borrowObject();
    md->copy( *this );
    return md;
}

void 
fluke::ModuleDownstream::copy( const ChromosomeElement &ce ) {
    Downstream::copy( ce );
    const ModuleDownstream *tp = 
        dynamic_cast< const ModuleDownstream * >( &ce );
    module_ = tp->module_;
}

void
fluke::ModuleDownstream::toPool() {
    ObjectCache< ModuleDownstream >::instance()->returnObject( this );
}

int 
fluke::ModuleDownstream::mutate() {
    return 0;
}

std::string 
fluke::ModuleDownstream::asXmlString() const {
    return "<dstream id=\"" + boost::lexical_cast< std::string >( tag_ ) + 
           "\" module=\"" + boost::lexical_cast< std::string >( module_ ) + 
           "\"/>\n";
}
