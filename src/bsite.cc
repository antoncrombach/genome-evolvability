//
// Implementation of the binding site object.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "bsite.hh"


fluke::ShortSeqManager* fluke::BindingSite::ssm_ = 0;
template<> fluke::ObjectCache< fluke::BindingSite >* 
fluke::ObjectCache< fluke::BindingSite >::instance_ = 0;


fluke::BindingSite::BindingSite() : ChromosomeElement(), tfbs_( -1 ) {}

fluke::BindingSite::BindingSite( label l ) : ChromosomeElement(), tfbs_( l ) {}

fluke::BindingSite::BindingSite( const BindingSite &bs ) 
    : ChromosomeElement( bs ) {
    BindingSite::copy( bs );
}

fluke::BindingSite::~BindingSite() {
    ssm_->freeShortSeq( tfbs_ );
}

fluke::ChromosomeElement*
fluke::BindingSite::clone() const {
    //return new BindingSite( *this );
    BindingSite *bs = ObjectCache< BindingSite >::instance()->borrowObject();
    bs->tfbs_ = tfbs_;
    // missing indirection does not matter
    ssm_->allocShortSeq( tfbs_ );
    return bs;
}

void 
fluke::BindingSite::copy( const ChromosomeElement &chr ) {
    const BindingSite *bs = dynamic_cast< const BindingSite * >( &chr );
    tfbs_ = bs->tfbs_;
    ssm_->allocShortSeq( tfbs_ );
}

void
fluke::BindingSite::toPool() {
    ObjectCache< BindingSite >::instance()->returnObject( this );
}

int 
fluke::BindingSite::mutate() {
    int result = 0;
    boost::tie( tfbs_, result ) = ssm_->mutateShortSeq( tfbs_ );
    return result;
}

void
fluke::BindingSite::shortSeqManager( ShortSeqManager *s ) 
{ ssm_ = s; }

fluke::ShortSeqManager*
fluke::BindingSite::shortSeqManager() 
{ return ssm_; }

std::string 
fluke::BindingSite::asXmlString() const {
    return "<bsite seq=\"" + ssm_->strShortSeq( tfbs_ ) + "\"/>";
}

std::string
fluke::BindingSite::asString() const {
    return ssm_->strShortSeq( tfbs_ );
}
