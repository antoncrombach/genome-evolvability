//
// Implementation of a genome agent with gene modules.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "module_agent.hh"
#include "population.hh"
#include "environment.hh"


float fluke::ModuleAgent::birth_rate_ = 0.0;
float fluke::ModuleAgent::death_rate_ = 0.0;

int fluke::ModuleAgent::max_dist_ = 0;
int fluke::ModuleAgent::penalty_genome_ = 0;
int fluke::ModuleAgent::penalty_tposons_ = 0;
double fluke::ModuleAgent::penalty_genome_rate_ = 0;
double fluke::ModuleAgent::penalty_tposons_rate_ = 0;
double fluke::ModuleAgent::dist_coeff_ = 0;

std::vector< fluke::Chromosome::tag_container >
fluke::ModuleAgent::module_tags_ = std::vector< Chromosome::tag_container >();

fluke::Chromosome::tag_container 
fluke::ModuleAgent::essential_tags_ = Chromosome::tag_container();


fluke::ModuleAgent::ModuleAgent( int tt, Genome *g ) 
    : Agent( tt ), mod_tags_now_(), ess_tags_now_() {
    genome_ = g;
    distance_ = 0;
    distance_parent_ = 0;
    size_parent_ = genome_->fullSize();
    inventorised_ = false;
}

fluke::ModuleAgent::ModuleAgent( const ModuleAgent &ag ) : Agent() {
    ModuleAgent::copy( ag );
}

fluke::ModuleAgent::~ModuleAgent() {
    delete genome_;
}

fluke::Agent* 
fluke::ModuleAgent::clone() const {
    ModuleAgent *aux = new ModuleAgent( *this );
    return aux;
}

void 
fluke::ModuleAgent::copy( const Agent &ag ) {
    const ModuleAgent *ma = dynamic_cast< const ModuleAgent * >( &ag );
    Agent::copy( *ma );
    distance_ = ma->distance_;
    distance_parent_ = ma->distance_parent_;
    size_parent_ = ma->size_parent_;
    genome_ = ma->genome_->clone();
    inventorised_ = false;
    // and copy the counts of genes
    std::copy( ma->ess_tags_now_.begin(), ma->ess_tags_now_.end(),
        std::back_inserter( ess_tags_now_ ) );
    for( const_module_iter ii = ma->mod_tags_now_.begin(); 
        ii != ma->mod_tags_now_.end(); ++ii ) {
        mod_tags_now_.push_back( std::vector< uint >( ii->begin(), ii->end() ));
    }
}

void
fluke::ModuleAgent::initialise() {
    // quick fix...
    if( essential_tags_.empty() ) referenceTags( *this );
    // Make sure the retroposons are counted
    genome_->nrRetroposons();
    // and that the length is cached
    genome_->fullSize();
}

void 
fluke::ModuleAgent::step( Population &pop ) {
    // only for testing purposes
    /*genome_->mutate();*/
    // only for testing purposes
    float rr = uniform();
    if( rr < death_rate_ ) {
        dying_ = true;
    }
}

fluke::Agent*
fluke::ModuleAgent::sibling() {
    Agent *result = 0;
    // update parent info
    distance_parent_ = distance_;
    size_parent_ = genome_->fullSize();
    // first perform sort-of mitosis
    genome_->duplicate();
    genome_->mutate();
    inventorised_ = false;
    Genome *sister_genome = genome_->split();
    // build sister agent
    ModuleAgent *sister = new ModuleAgent( type_, sister_genome );
    std::copy( ess_tags_now_.begin(), ess_tags_now_.end(),
        std::back_inserter( sister->ess_tags_now_ ) );
    for( module_iter ii = mod_tags_now_.begin(); 
        ii != mod_tags_now_.end(); ++ii ) {
        sister->mod_tags_now_.push_back( 
            std::vector< uint >( ii->begin(), ii->end() ) );
    }
    // update parent info
    sister->distance_parent_ = distance_parent_;
    sister->size_parent_ = size_parent_;
    // return
    result = sister;    
    return result;
}

double
fluke::ModuleAgent::score() const {
    // number of genes still to be copied is already in distance_
    double result = static_cast< double >( distance_ );
    // now see if we have any penalties to add
    int aux = genome_->fullSize() - penalty_genome_;
    int bux = genome_->nrRetroposons() - penalty_tposons_;
    if( aux > 0 ) {
        result += aux * penalty_genome_rate_;
    }
    if( bux > 0 ) {
        result += bux * penalty_tposons_rate_;
    }
    // reviewer 2 centromere check
    if( !genome_->oneCentromere() ) {
        result += max_dist_;
    }
    // we've got the raw score now and calculate the fitness score
    if( result < max_dist_ ) {
        return 1.0 - dist_coeff_ * result;
    } else {
        return 0.0;
    }
}

void
fluke::ModuleAgent::evaluate( const Environment &env ) {
    if( !inventorised_ ) { 
        countGenes();
    }
    distance_ = essentialsScore( env ) + modulesScore( env );
}

double
fluke::ModuleAgent::penalty( double ss ) const {
    return ss;
}

void
fluke::ModuleAgent::countGenes() {
    // inventorise for essential genes and all modules
    countEssentialGenes();
    countModuleGenes();
    inventorised_ = true;
}

void
fluke::ModuleAgent::countEssentialGenes() {
    // get the essential genes
    Chromosome::tag_container aux = genome_->essentialTags();
    std::sort( aux.begin(), aux.end() );
    // loop through two sorted vectors
    ess_tags_now_.clear();
    std::fill_n( std::back_inserter( ess_tags_now_ ), 
        essential_tags_.size(), 0 );
    uint i = 0;
    Chromosome::tag_container::iterator j = aux.begin();
    while( i != essential_tags_.size() ) {
        if( j != aux.end() ) {
            if( *j == essential_tags_[ i ] ) {
                ++ess_tags_now_[ i ];
                ++j;
            } else {
                ++i;
            }
        } else {
            i = essential_tags_.size();
        }
    }    
}

void
fluke::ModuleAgent::countModuleGenes() {
    int jj = 0;
    module_iter ii = module_tags_.begin();
    if( mod_tags_now_.empty() ) {
        std::fill_n( std::back_inserter( mod_tags_now_ ), 
            module_tags_.size(), Chromosome::tag_container() );
    }
    while( ii != module_tags_.end() ) {
        // see which genes we have
        Chromosome::tag_container aux = genome_->moduleTags( jj );
        std::sort( aux.begin(), aux.end() );
        // Loop through tags
        mod_tags_now_[ jj ].clear();
        std::fill_n( std::back_inserter( mod_tags_now_[ jj ] ), ii->size(), 0 );
        uint i = 0;
        Chromosome::tag_container::iterator j = aux.begin();
        while( i != ii->size() ) {
            if( j != aux.end() ) {
                if( *j == ( *ii )[ i ] ) {
                    ++mod_tags_now_[ jj ][ i ];
                    ++j;
                } else {
                    ++i;
                }
            } else {
                i = ii->size();
            }
        }
        ++ii;
        ++jj;
    }
}

int
fluke::ModuleAgent::essentialsScore( const Environment &env ) {
    // expected nr of copies is one?
    int result = 0;
    uint i = 0;
    while( i != ess_tags_now_.size() ) {
        // minus one, because its the distance from having one copy
        int aux = ess_tags_now_[ i ] - 1;
        if( aux >= 0 && aux < max_dist_ ) {
            result += aux;
            ++i;
        } else {
            result = max_dist_;
            i = ess_tags_now_.size();
        }
    }
    return result;
}

int
fluke::ModuleAgent::modulesScore( const Environment &env ) {
    if( mod_tags_now_.empty() ) return max_dist_;
    // expected nr of copies depends on environment
    int result = 0;
    uint i = 0;
    while( i != mod_tags_now_.size() ) {
        uint j = 0;
        while( j != mod_tags_now_[ i ].size() ) {
            if( mod_tags_now_[ i ][ j ] != 0 ) {
                int aux = abs( mod_tags_now_[ i ][ j ] 
                    - env.expectedCopies( i ) );
                if( aux >= 0 && aux < max_dist_ ) {
                    result += aux;
                }
                ++j;
            } else {
                j = mod_tags_now_[ i ].size();
                i = mod_tags_now_.size() - 1;
                result = max_dist_;
            }
        }
        ++i;
    }
    return result;
}

void 
fluke::ModuleAgent::write( std::ostream &os ) const {
    os << "<agent type=\"" << type_ << "\" birth=\"" << me_.time;
    os << "\" x=\"" << me_.x << "\" y=\"" << me_.y << "\" i=\"" << me_.i;
    os << "\">\n<class>ModuleAgent</class>\n";
    os << "<score>" << score() << "</score>\n";
    os << "<parent time=\"" << ancestor_.time 
       << "\" x=\"" << ancestor_.x << "\" y=\"" << ancestor_.y
       << "\" i=\"" << ancestor_.i << "\"/>\n";
    os << *genome_; 
    os << "<mods>";
    for( const_module_iter i = mod_tags_now_.begin(); 
        i != mod_tags_now_.end(); ++i ) {
        os <<"<mod>";
        std::copy( i->begin(), i->end(), 
            std::ostream_iterator< uint >( os, " " ) );
        os << "</mod>\n";
    }
    os << "</mods>\n";

    os << "<ess>";
    std::copy( ess_tags_now_.begin(), ess_tags_now_.end(),
        std::ostream_iterator< uint >( os, " " ) );
    os << "</ess>\n";
    //os << "<debug retro =\"" << genome_->nrRetroposons() << "\"/>\n";
    os << "<genodist>" << distance_ << "</genodist>\n";
    os << "</agent>\n";
}

void
fluke::ModuleAgent::birthRate( float f )
{ birth_rate_ = f; }

float
fluke::ModuleAgent::birthRate()
{ return birth_rate_; }

void
fluke::ModuleAgent::deathRate( float f )
{ death_rate_ = f; }

float
fluke::ModuleAgent::deathRate()
{ return death_rate_; }

void
fluke::ModuleAgent::referenceTags( const ModuleAgent &ma ) {
    // check which genes we have and create a reference lookup
    std::vector< uint > aux = ma.genome_->essentialTags();
    std::sort( aux.begin(), aux.end() );
    Chromosome::tag_iter last = std::unique( aux.begin(), aux.end() );
    std::copy( aux.begin(), last, std::back_inserter( essential_tags_ ) );
    aux.clear();
    // and now the modules
    int i = 0;
    aux = ma.genome_->moduleTags( i );
    // empty module implies no more modules present
    while( !aux.empty() ) {
        std::sort( aux.begin(), aux.end() );
        last = std::unique( aux.begin(), aux.end() );
        module_tags_.push_back( std::vector< uint >( aux.begin(), last ) );
        aux.clear();
        aux = ma.genome_->moduleTags( ++i );
    }
}

void
fluke::ModuleAgent::maxDistance( int f ) 
{ max_dist_ = f; dist_coeff_ = 1.0 / max_dist_; }

int
fluke::ModuleAgent::maxDistance() 
{ return max_dist_; }

void
fluke::ModuleAgent::maxGenomeSize( int f )
{ penalty_genome_ = f; }

void
fluke::ModuleAgent::maxTposons( int f )
{ penalty_tposons_ = f; }

void
fluke::ModuleAgent::genomePenaltyRate( double f )
{ penalty_genome_rate_ = f; }

void
fluke::ModuleAgent::retroposonPenaltyRate( double f )
{ penalty_tposons_rate_ = f; }
