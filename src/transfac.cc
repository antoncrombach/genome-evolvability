//
// Implementation of a transcription factor downstream region.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "transfac.hh"


fluke::ShortSeqManager* fluke::TranscriptionFactor::ssm_ = 0;
template<> fluke::ObjectCache< fluke::TranscriptionFactor >* 
fluke::ObjectCache< fluke::TranscriptionFactor >::instance_ = 0;


fluke::TranscriptionFactor::TranscriptionFactor( 
    const TranscriptionFactor &tp ) {
    TranscriptionFactor::copy( tp );
}

fluke::TranscriptionFactor::~TranscriptionFactor() {
    ssm_->freeShortSeq( tf_ );
}

fluke::ChromosomeElement* 
fluke::TranscriptionFactor::clone() const {
    //return new TranscriptionFactor( *this );
    TranscriptionFactor *tp = 
        ObjectCache< TranscriptionFactor >::instance()->borrowObject();
    tp->tf_ = tf_;
    ssm_->allocShortSeq( tf_ );
    return tp;
}

void 
fluke::TranscriptionFactor::copy( const ChromosomeElement &ce ) {
    Downstream::copy( ce );
    const TranscriptionFactor *tp = 
        dynamic_cast< const TranscriptionFactor * >( &ce );
    tf_ = tp->tf_;
    ssm_->allocShortSeq( tf_ );
}

void
fluke::TranscriptionFactor::toPool() {
    ObjectCache< TranscriptionFactor >::instance()->returnObject( this );
}

int 
fluke::TranscriptionFactor::mutate() {
    int result = 0;
    boost::tie( tf_, result ) = ssm_->mutateShortSeq( tf_ );
    return result;
}

void
fluke::TranscriptionFactor::shortSeqManager( ShortSeqManager *s ) 
{ ssm_ = s; }

std::string 
fluke::TranscriptionFactor::asXmlString() const {
    return "<transfac id=\"" + boost::lexical_cast< std::string >( tag_ ) + 
           "\">" + ssm_->strShortSeq( tf_ ) + "</transfac>\n";
}
