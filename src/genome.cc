//
// Genome class implementation.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "genome.hh"
#include "chromosome.hh"
#include "chromelement.hh"

fluke::Genome::Genome() : chromos_( new std::list< Chromosome* >() ),
    nr_dsbs_parent_(), nr_mutations_() {}

// note: using magic number
fluke::Genome::Genome( std::list< Chromosome* > *c ) 
: chromos_( c ), nr_dsbs_parent_( 2 * c->size(), 0 ), nr_mutations_( 6, 0 ) {
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).parent( this );
    }
}

fluke::Genome::Genome( const Genome &g ) : nr_dsbs_parent_(), nr_mutations_() {
    Genome::copy( g );
}

fluke::Genome::~Genome() {
    clear();
    delete chromos_;
}

fluke::Genome*
fluke::Genome::clone() const {
    return new Genome( *this );
}

void
fluke::Genome::copy( const Genome &g ) {
    chromos_ = new std::list< Chromosome* >();
    for( chromos_iter i = g.chromos_->begin(); i != g.chromos_->end(); ++i ) {
        //chromos_->push_back( new Chromosome( **i ) );
        chromos_->push_back( ( **i ).clone() );
        chromos_->back()->parent( this );
    }
    std::copy( g.nr_dsbs_parent_.begin(), g.nr_dsbs_parent_.end(),
        std::back_inserter( nr_dsbs_parent_ ) );
    std::copy( g.nr_mutations_.begin(), g.nr_mutations_.end(),
        std::back_inserter( nr_mutations_ ) );
}

void
fluke::Genome::duplicate() {
    std::list< Chromosome* > aux;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        //aux.push_back( new Chromosome( **i ) );
        aux.push_back( ( **i ).clone() );
    }
    chromos_->splice( chromos_->end(), aux );
}

fluke::Genome*
fluke::Genome::split() {
    // note: assuming only two chromosomes
    std::list< Chromosome* > *aux = new std::list< Chromosome* >();
    aux->splice( aux->begin(), *chromos_, boost::next( chromos_->begin() ) );
    Genome *result = new Genome( aux );
    std::copy( nr_dsbs_parent_.begin(), nr_dsbs_parent_.end(),
        result->nr_dsbs_parent_.begin() );
    std::copy( nr_mutations_.begin(), nr_mutations_.end(),
        result->nr_mutations_.begin() );
    return result;
}

int 
fluke::Genome::mutate() {
    int result = 0;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        ( **i ).mutate();
    }
    // Gather a few counts, first clear
    nr_dsbs_parent_.clear();
    std::fill( nr_mutations_.begin(), nr_mutations_.end(), 0 );
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        // genes
        std::vector< uint > aux = ( **i ).nrMutations();
        // might be empty if a cached version is taken...
        if( ! aux.empty() ) {
            std::transform( nr_mutations_.begin(), nr_mutations_.end(),
                aux.begin(), nr_mutations_.begin(), std::plus< uint >() );
        }
        // dsbs
        nr_dsbs_parent_.push_back( ( **i ).nrDoubleStrandBreaks() );
    }
    
    // get all segments
    std::list< Chromosome* > recombined;
    std::vector< Chromosome* > heads, middles, tails;
    heads.reserve( chromos_->size() );
    tails.reserve( chromos_->size() );
    chromos_iter i = chromos_->begin();
    while( i != chromos_->end() ) {
        // get the segments
        std::list< Chromosome* > aux = ( **i ).segments();
        // do something depending on # segments
        if( aux.size() == 1 ) {
            recombined.push_back( aux.front() );
        } else if( aux.size() == 2 ) {
            heads.push_back( aux.front() );
            tails.push_back( aux.back() );
        } else {
            middles.reserve( middles.size() + aux.size() - 2 );
            heads.push_back( aux.front() );
            tails.push_back( aux.back() );
            // and all the middle parts
            std::copy( boost::next( aux.begin() ), boost::prior( aux.end() ),
                std::back_inserter( middles ) );
        }
        aux.clear();
        ++i;
    }

    // try to reuse these...
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        //if( ( **i ).empty() ) delete *i;
        if( ( **i ).empty() ) ( **i ).toPool();
    }

    // randomly assign middle segments to heads
    std::random_shuffle( middles.begin(), middles.end(), rand_range< int > );
    for( std::vector< Chromosome* >::iterator i = middles.begin(); 
        i != middles.end(); ++i ) {
        Chromosome *aux = *( random_element( heads.begin(),
             heads.end(), rand_range< int > ) );
        aux->append( *i );
    }
    //smart_erase( middles, middles.begin(), middles.end() );
    smart_return( middles, middles.begin(), middles.end() );

    // randomly assign the last part of the chromosomes
    std::random_shuffle( tails.begin(), tails.end(), rand_range< int > );
    for( std::vector< Chromosome* >::iterator i = heads.begin(); 
        i != heads.end(); ++i ) {
        ( **i ).append( tails.back() );
        delete tails.back();
        tails.pop_back();
    }

    // make them recache their info
    for( std::vector< Chromosome* >::iterator i = heads.begin(); 
        i != heads.end(); ++i ) {
        ( **i ).recache();
    }

    // and copy the new ones 
    chromos_->clear();
    chromos_->splice( chromos_->begin(), recombined );
    std::copy( heads.begin(), heads.end(), std::back_inserter(
            *chromos_ ) );
    heads.clear();

    return result; 
}

fluke::Chromosome::tag_container
fluke::Genome::essentialTags() const {
    Chromosome::tag_container result;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        Chromosome::tag_container aux = ( **i ).essentialTags();
        // valid because we have only one chromosome...
        result.insert( result.end(), aux.begin(), aux.end() );
    }
    return result;
}

fluke::Chromosome::tag_container
fluke::Genome::moduleTags( int module ) const {
    Chromosome::tag_container result;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        Chromosome::tag_container aux = ( **i ).moduleTags( module );
        // valid because we have only one chromosome...
        result.insert( result.end(), aux.begin(), aux.end() );
    }
    return result;
}

boost::tuple< fluke::Chromosome*, 
    std::list< fluke::ChromosomeElement* >::iterator >
fluke::Genome::randElement() {
    Chromosome *bux = *( random_element( chromos_->begin(), chromos_->end(),
        rand_range< int > ) );
    return bux->randChromosomeElement();
}

void
fluke::Genome::write( std::ostream &os ) const {
    os << "<genome>\n";
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        os << **i;
    }
    os << "</genome>\n";
}

int
fluke::Genome::fullSize() const {
    int result = 0;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).size();
    }
    return result;
}

void
fluke::Genome::clear() {
    //smart_erase( *chromos_, chromos_->begin(), chromos_->end() );
    smart_return( *chromos_, chromos_->begin(), chromos_->end() );
}

int
fluke::Genome::nrRetroposons() const {
    int result = 0;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrRetroposons();
    }
    return result;
}

int
fluke::Genome::nrRepeats() const {
    int result = 0;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrRepeats();
    }
    return result;
}

int
fluke::Genome::nrDoubleStrandBreaks() const {
    int result = 0;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result += ( **i ).nrDoubleStrandBreaks();
    }
    return result;
}

bool
fluke::Genome::oneCentromere() const {
    bool result = true;
    for( chromos_iter i = chromos_->begin(); i != chromos_->end(); ++i ) {
        result = result && ( **i ).oneCentromere();
    }
    return result;
}
