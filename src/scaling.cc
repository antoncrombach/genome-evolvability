//
// Implementation of some scaling schemes.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "scaling.hh"

void 
fluke::NoScaling::scale( std::vector< double > &scores ) {
    // skip
}

void 
fluke::LinearScaling::scale( std::vector< double > &scores ) {
    // pre: scores are in range [0..++)
    //
    // note: as in the previous Python implementation, I patched the algorithm 
    // with the feature that if the max is zero, all the scores are set to 
    // a base fitness (otherwise all the scores would be zero and the 
    // selection would always pick the not-do-anything option).
    // get the max
    double mx = *( std::max_element( scores.begin(), scores.end() ) );
    // if the max is almost zero, they're all zero. So fill with base score
    if( close_to( mx, static_cast< double >( 0.0 ) ) ) {
        std::fill( scores.begin(), scores.end(), base_score_ );
    } else {
        // note: if all scores are the same, the scaling will result in all
        // scores being one
        std::transform( scores.begin(), scores.end(), 
            scores.begin(), std::bind2nd( std::divides< double >(), mx ) );
    }
}

void 
fluke::PowerScaling::scale( std::vector< double > &scores ) {
    // pre: scores are in range [0..1)
    //
    // note: as in the previous Python implementation, I patched the algorithm 
    // with the feature that if the max is zero, all the scores are set to 
    // a base fitness (otherwise all the scores would be zero and the 
    // selection would always pick the not-do-anything option).
    // get the max
    double mx = *( std::max_element( scores.begin(), scores.end() ) );
    // if the max is almost zero, they're all zero. So fill with base score
    if( close_to( mx, static_cast< double >( 0.0 ) ) ) {
        std::fill( scores.begin(), scores.end(), base_score_ );
    } else {
        std::transform( scores.begin(), scores.end(), 
            scores.begin(), std::bind2nd( std::ptr_fun( pow ), power_ ) );
    }
}

