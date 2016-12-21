//
// Implementation of selection schemes.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "selection.hh"

fluke::Agent* 
fluke::RandSelection::select(
        std::vector< Agent* > &va, std::vector< double > &sc ) {
    // note: ignoring scores
    return va[ static_cast< int >( rand_range( va.size() ) ) ];
}

fluke::Agent* 
fluke::ProbalisticSelection::select( 
        std::vector< Agent* > &va, std::vector< double > &sc ) {
    // note: the "no-selection" case may be a 0-pointer in 'va'
    std::partial_sum( sc.begin(), sc.end(), sc.begin() );
    std::vector< double >::iterator aux = 
        std::upper_bound( sc.begin(), sc.end(), rand_range( sc.back() ) );
    int i = std::distance( sc.begin(), aux );
    return va[ i ];
}

