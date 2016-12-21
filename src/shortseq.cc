//
// Implementation of the short sequences flyweight manager.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "shortseq.hh"

using fluke::generator;

int fluke::ShortSeqManager::max_hamming_ = 0;

fluke::ShortSeqManager::ShortSeqManager( std::string alphabet,
    int length, double mutation_rate ) :
    rand_length_( generator, boost::uniform_int<>( 0, length - 1 ) ),
    rand_repository_( generator, boost::uniform_int<>( 0, 0 ) ),
    rand_alphabet_( generator, boost::uniform_int<>( 0, 0 ) ) {
    // initialise all the attributes
    sort( alphabet.begin(), alphabet.end() );
    unique_copy( alphabet.begin(), alphabet.end(), back_inserter( alphabet_ ));
    length_ = length;
    mutation_rate_ = 1.0 - pow( 1.0 - mutation_rate, length_ );
    
    // initialise some random nr fluke::generators correctly
    rand_repository_ = randrange_gen_type( fluke::generator, 
            boost::uniform_int<>( 0, static_cast< int >( 
                    pow( alphabet_.length(), length_ ) ) - 1 ) ),
    rand_alphabet_ = randrange_gen_type( fluke::generator, 
            boost::uniform_int<>( 1, alphabet_.length() - 1 ) ),

    // create all shortseqs
    generateShortSeqs();
}


fluke::ShortSeqManager::~ShortSeqManager() {
    delete short_seqs_;
    delete references_;
}

std::string
fluke::ShortSeqManager::strShortSeq( label lbl ) const {
    // return the std::string of this integer encoding
#ifdef DEBUG
    return short_seqs_->at( lbl );
#else
    return ( *short_seqs_ )[ lbl ];
#endif
}

void
fluke::ShortSeqManager::allocShortSeq( label lbl ) {
    // allocate another instance of the int encoded sequence
#ifdef DEBUG
    references_->at( lbl )++;
#else
    ( *references_ )[ lbl ]++;
#endif
}

void
fluke::ShortSeqManager::freeShortSeq( label lbl ) {
    // notify that int encoded sequence is no longer used
#ifdef DEBUG
    references_->at( lbl )--;
#else
    ( *references_ )[ lbl ]--;
#endif
}

fluke::label
fluke::ShortSeqManager::allocRandShortSeq() {
    // randomly get a short sequence from the repository, returned as label
    label rr = rand_repository_();
    allocShortSeq( rr );
    return rr;
}

boost::tuple< fluke::label, int >
fluke::ShortSeqManager::mutateShortSeq( int lbl ) {
    // apply point mutations on this sequence and return the new one
    if( close_to( mutation_rate_, 0.0 ) ) return boost::make_tuple( lbl, 0 );
#ifdef DEBUG
    std::string tmp = short_seqs_->at( lbl );
#else
    std::string tmp = ( *short_seqs_ )[ lbl ];
#endif
    int aux = 0;
    label result = lbl;
    if( uniform() < mutation_rate_ ) {
        freeShortSeq( lbl );
        int i = rand_length_();
        std::string a = alphabet_;
        int j = a.find_first_of( tmp[ i ] );
        
        // get random (but diff) nucleotide
        std::swap( a[ 0 ], a[ j ] );
        tmp[ i ] = a[ rand_alphabet_() ];
        // binary search for short seq's label and calculate index
        result = std::distance( short_seqs_->begin(),
            lower_bound( short_seqs_->begin(), short_seqs_->end(), tmp ) );
        allocShortSeq( result );
        aux = 1;
    }
    return boost::make_tuple( result, aux ); 
}

bool
fluke::ShortSeqManager::similar( label a, label b ) {
    std::string aux = ( *short_seqs_ )[ a ];
    std::string bux = ( *short_seqs_ )[ b ];

    int count = 0;
    std::string::iterator i = aux.begin();
    std::string::iterator j = bux.begin();
    while( i != aux.end() ) {
        if( *i == *j ) {
            ++count;
        }
        ++i;
        ++j;
    }
    return length_ - count <= max_hamming_;
}

const std::vector< int >& 
fluke::ShortSeqManager::referenceCounts() const {
    // return the reference counts
    return *references_;
}

void 
fluke::ShortSeqManager::generateShortSeqs() {
    // generate all the sequences, given the alphabet and the length
    int n = static_cast< int >( pow( alphabet_.length(), length_ ) );
    std::string i( length_, alphabet_[ 0 ] );

    short_seqs_ = new std::vector< std::string >( n, i );
    references_ = new std::vector< int >( n, 0 );
   
    // this is similar to counting
    for( std::vector< std::string >::iterator j = short_seqs_->begin() + 1; 
            j != short_seqs_->end(); j++ ) {
        *j = *( j - 1 );

        int tmp = j->find_last_not_of( 
            alphabet_.substr( alphabet_.length() - 1, alphabet_.length() ) );
        if( tmp != length_ - 1 ) {
            j->replace( tmp + 1, length_ - tmp - 1, 
                length_ - tmp - 1, alphabet_[ 0 ] );
        }
        ( *j )[ tmp ] = alphabet_[ alphabet_.find( ( *j )[ tmp ] ) + 1 ]; 
    }
}

void
fluke::ShortSeqManager::maxDistance( int f )
{ max_hamming_ = f; }

int
fluke::ShortSeqManager::maxDistance()
{ return max_hamming_; }

