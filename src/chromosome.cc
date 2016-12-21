//
// Chromosome implementation...
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "pool.hh"
#include "chromosome.hh"


template<> fluke::ObjectCache< fluke::Chromosome >* 
fluke::ObjectCache< fluke::Chromosome >::instance_ = 0;
fluke::MutateRates *fluke::Chromosome::rate_mutator_ = 0;

// note: using magic number
fluke::Chromosome::Chromosome() 
    : parent_( 0 ), chro_( new std::list< ChromosomeElement* >() ),
      mut_events_( 6, 0 ), nr_retroposons_( 0 ), nr_ltr_( 0 ), len_( 0 ),
      update_retro_( true ), update_ltr_( true ), update_len_( true ) {}

// note: using magic number
fluke::Chromosome::Chromosome( Genome *g, std::list< ChromosomeElement* > *ll ) 
    : parent_( g ), chro_( ll ), mut_events_( 6, 0 ), nr_retroposons_( 0 ),
      nr_ltr_( 0 ), len_( ll->size() ), update_retro_( true ),
      update_ltr_( true ), update_len_( true ) {}

fluke::Chromosome::Chromosome( const Chromosome &c ) 
    : chro_( new std::list< ChromosomeElement* >() ), mut_events_( 6, 0 ) {
    copy( c );
}

fluke::Chromosome::~Chromosome() {
    smart_erase( *chro_, chro_->begin(), chro_->end() );
    delete chro_;
}

void
fluke::Chromosome::copy( const Chromosome &c ) {
    parent_ = c.parent_;
    for( ce_iter i = c.chro_->begin(); i != c.chro_->end(); ++i ) {
        ChromosomeElement *aux = ( **i ).clone();
        chro_->push_back( aux );
    }
    std::copy( c.mut_events_.begin(), c.mut_events_.end(), mut_events_.begin());
    nr_retroposons_ = c.nr_retroposons_;
    nr_ltr_ = c.nr_ltr_;
    len_ = c.len_;
    update_retro_ = c.update_retro_;
    update_ltr_ = c.update_ltr_;
    update_len_ = c.update_len_;
    copyRates( c, *this );
}

fluke::Chromosome*
fluke::Chromosome::clone() const {
    Chromosome *cc = ObjectCache< Chromosome >::instance()->borrowObject();
    cc->copy( *this );
    return cc;
}

void
fluke::Chromosome::toPool() {
    // empty everything
    smart_return( *chro_, chro_->begin(), chro_->end() );
    std::fill_n( mut_events_.begin(), 6, 0 );
    nr_retroposons_ = 0;
    nr_ltr_ = 0;
    len_ = 0;    
    update_retro_ = true;
    update_ltr_ = true;
    update_len_ = true;
    ObjectCache< Chromosome >::instance()->returnObject( this );
}

int
fluke::Chromosome::mutate() {
    // make sure everything is in correct state
    reset();
    // first insert new retroposons
    newRetrotransposon();
    // now loop
    ce_iter i = chro_->begin();
    while( i != chro_->end() ) {
        // note: it is a bit tricky, but incrementing i is done in the methods
        if( ( **i ).isActive() ) {
            /* if( IsBindingSite()( *i ) ) {
                result += ( **i ).mutate();
                i = bsiteMutate( i );
            } else */ if( IsTrueDstream()( *i ) ) {
                // note: downstreams do not have mutations internally
                /* result += ( **i ).mutate(); */
                i = geneMutate( i );
            } else if( IsRepeat()( *i ) ) {
                i = repeatMutate( i );
            } else if( IsRetroposon()( *i ) ) {
                i = retroposonMutate( i );
            } else if( IsCentromere()( *i ) ) {
                ++i;
            }
        } else {
            ++i;
        }
    }
    // and now mutate rates
    if( !close_to( mut_rate_, 0.0 ) ) rate_mutator_->mutate( *this );
    return 0;
}

void
fluke::Chromosome::newRetrotransposon() {
    if( close_to( new_tp_rate_, 0.0 ) ) {
        return;
    }
    // FIX using magic number!!!
    // assuming small rate
    int aux = size();
    if( uniform() < 1.0 - pow( 1.0 - new_tp_rate_, aux ) ) {
        // find a spot
        Chromosome *bux;
        ce_iter cux;
        boost::tie( bux, cux ) = randGenomeElement();
        // insert
        bux->chro_->insert( cux, new Repeat() );
        bux->chro_->insert( cux, new Retroposon( rand_range( 100 ) + 100 ) );
        bux->chro_->insert( cux, new Repeat() );
        ++nr_retroposons_;
        nr_ltr_ += 2;
        len_ += 3;
    }
}

fluke::Chromosome::ce_iter
fluke::Chromosome::repeatMutate( ce_iter i ) {
    if( close_to( dsb_recombination_ + rm_ltr_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    Repeat *aux = dynamic_cast< Repeat* >( *i );
    // just in case (patch)
    aux->repairDSB();
    double uu = uniform();
    if( uu < dsb_recombination_ ) {
        // needs to be repaired @ genome level
        aux->induceDSB();
        ( **i ).inactivate();
        ++i;
        ++mut_events_[ DSB ];
        //std::cout << "dsb" << std::endl;
    } else if( uu < dsb_recombination_ + rm_ltr_rate_ ) {
        // first check if flanking a retroposon
        if( i == chro_->begin() ) {
            if( IsRetroposon()( *( boost::next( i ) ) ) ) {
                ++i;
            } else {
                i = smart_return( *chro_, i );
                ++mut_events_[ RM_LTR ];
                //std::cout << "ltr" << std::endl;
                --nr_ltr_;
                --len_;
            }
        } else if( boost::next( i ) == chro_->end() ) {
            if( IsRetroposon()( *( boost::prior( i ) ) ) ) {
                ++i;
            } else {
                i = smart_return( *chro_, i );
                ++mut_events_[ RM_LTR ];
                //std::cout << "ltr" << std::endl;
                --nr_ltr_;
                --len_;
            }
        } else {
            if( IsRetroposon()( *( boost::prior( i ) ) ) ||
                IsRetroposon()( *( boost::next( i ) ) ) ) {
                ++i;
            } else {
                i = smart_return( *chro_, i );
                ++mut_events_[ RM_LTR ];
                //std::cout << "ltr" << std::endl;
                --nr_ltr_;
                --len_;
            }
        }
    } else {
        ++i;
    }
    return i;
}

fluke::Chromosome::ce_iter
fluke::Chromosome::bsiteMutate( ce_iter i ) {
    /*// FIX: not keeping track of length here!
    // pre: i is active and binding site
    if( close_to( nw_bs_rate_ + cp_bs_rate_ + rm_bs_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double rr = uniform();
    if( rr < nw_bs_rate_ ) {
        // insert a new bsite somewhere in the genome
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeElement();
        label cux = BindingSite::shortSeqManager()->allocRandShortSeq();
        ce_iter j = aux->insert( bux, new BindingSite( cux ) );
        ( **j ).inactivate();
        ++i;
    } else if( rr < nw_bs_rate_ + cp_bs_rate_ ) {
        // insert a copy of the current bsite somewhere
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeElement();
        ce_iter j = aux->insert( bux, ( **i ).clone() );
        ( **j ).inactivate();
        ++i;
    } else if( rr < nw_bs_rate_ + cp_bs_rate_ + rm_bs_rate_ ) {
        // delete the current bsite
        i = smart_return( *chro_, i );
        --len_;
    } else {
        ++i;
    }*/
    return i;
}

fluke::Chromosome::ce_iter
fluke::Chromosome::geneMutate( ce_iter i ) {
    // pre: i is active and ( ordinary or module downstream )
    if( close_to( cp_gene_rate_ + rm_gene_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double uu = uniform();
    if( uu < cp_gene_rate_ ) {
        // insert a copy of the current gene somewhere in the genome
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeElement();
        // find the gene
        ce_iter ll = upstreamSelect( i );
        ce_iter rr = boost::next( i );
        std::list< ChromosomeElement* > cux;
        uint tt = 0;
        for( ce_iter j = ll; j != rr; ++j ) {
            ChromosomeElement *dux = ( **j ).clone();
            dux->inactivate();
            cux.push_back( dux );
            ++tt;
        }
        aux->splice( bux, cux, tt, 0 );
        ++mut_events_[ CP_G ];
        //std::cout << "cpg" << std::endl;
        ++i;
    } else if( uu < cp_gene_rate_ + rm_gene_rate_ ) {
        // delete the current gene
        ce_iter ll = upstreamSelect( i );
        ce_iter rr = boost::next( i );
        len_ -= std::distance( ll, rr );
        i = smart_return( *chro_, ll, rr );
        ++mut_events_[ RM_G ];
        //std::cout << "rmg" << std::endl;
    } else {
        ++i;
    }
    return i;
}

fluke::Chromosome::ce_iter
fluke::Chromosome::retroposonMutate( ce_iter i ) {
    // pre: i is active and retroposon
    if( close_to( cp_tp_rate_ + rm_tp_rate_, 0.0 ) ) {
        return boost::next( i );
    }
    double uu = uniform();
    if( uu < cp_tp_rate_ ) {
        // insert a copy of the retroposon and LTRs somewhere
        Chromosome *aux;
        ce_iter bux;
        boost::tie( aux, bux ) = randGenomeElement();
        // get the repeats
        ce_iter ll = boost::prior( i );
        ce_iter rr = boost::next( i, 2 );
        std::list< ChromosomeElement* > cux;
        for( ce_iter j = ll; j != rr; ++j ) {
            ChromosomeElement *dux = ( **j ).clone();
            dux->inactivate();
            cux.push_back( dux );
        }
        aux->splice( bux, cux, 3, 1 );
        ++mut_events_[ CP_RP ];
        //std::cout << "cprp" << std::endl;
        ++i;
    } else if( uu < cp_tp_rate_ + rm_tp_rate_ ) {
        // reciprocal recombination, one LTR stays
        ce_iter ll = boost::prior( i );
        ce_iter rr = boost::next( i );
        i = smart_return( *chro_, ll, rr );
        // inactivate leftover repeat
        ( **i ).inactivate();
        ++mut_events_[ RM_RP ];
        //std::cout << "rmrp" << std::endl;
        --nr_retroposons_;
        --nr_ltr_;
        len_ -= 2;
    } else {
        ++i;
    }
    return i;
}

std::list< fluke::Chromosome* >
fluke::Chromosome::segments() {
    std::list< Chromosome* > result;
    if( mut_events_[ DSB ] == 0 ) {
        // creating an alias...
        result.push_back( this );
    } else {
        // creation of new chromosomes automatically sets their update flag
        // get the bits and pieces
        bool found = false;
        ce_iter jj( chro_->begin() );
        ce_iter ii( jj );
        // and loop...
        while( jj != chro_->end() ) {
            Repeat *aux = dynamic_cast< Repeat* >( *jj );
            if( aux ) {
                if( aux->hasDSB() ) {
                    aux->repairDSB();
                    found = true;
                }
            }
            ++jj;
            // if found, put the segment in the vector
            if( found ) {
                /*
                std::list< ChromosomeElement* > *bux = 
                    new std::list< ChromosomeElement* >();
                bux->splice( bux->end(), *chro_, ii, jj );
                Chromosome *cux = new Chromosome( parent_, bux );
                */
                Chromosome *cux = 
                    ObjectCache< Chromosome >::instance()->borrowObject();
                cux->parent_ = parent_;
                cux->chro_->splice( cux->chro_->end(), *chro_, ii, jj );
                /*
                std::fill_n( std::back_inserter( cux->mut_events_ ), 6, 0 );
                */
                copyRates( *this, *cux );
                result.push_back( cux );
                ii = jj;
                found = false;
            }
        }
        // and the last one
        /*
        std::list< ChromosomeElement* > *bux = 
            new std::list< ChromosomeElement* >();
        bux->splice( bux->end(), *chro_, ii, jj );
        Chromosome *cux = new Chromosome( parent_, bux );
        */
        Chromosome *cux = 
            ObjectCache< Chromosome >::instance()->borrowObject();
        cux->parent_ = parent_;
        cux->chro_->splice( cux->chro_->end(), *chro_, ii, jj );
        /*
        std::fill_n( std::back_inserter( cux->mut_events_ ), 6, 0 );
        */
        copyRates( *this, *cux );
        result.push_back( cux );
    }
    // an empty chromosome is the leftover, if there were segments
    return result;
}

void
fluke::Chromosome::append( Chromosome *chr ) {
    // not copying parent_
    chro_->splice( chro_->end(), *( chr->chro_ ) );
    //++mut_events_[ DSB ];
}

fluke::Chromosome::ce_iter
fluke::Chromosome::insert( ce_iter i, ChromosomeElement *ce ) {
    ++len_;
    return chro_->insert( i, ce );
}

void
fluke::Chromosome::splice( ce_iter i, 
    std::list< ChromosomeElement* > ces, uint ll, uint rr ) {
    // overloading list method
    len_ += ll;
    nr_retroposons_ += rr;
    nr_ltr_ += 2 * rr;
    chro_->splice( i, ces );
}

void
fluke::Chromosome::reset() {
    // pre: chro_ is initialised
    for( ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        ( **i ).activate();
    }
    // making sure all entries exist and are zero
    std::fill_n( mut_events_.begin(), 6, 0 );
}

void
fluke::Chromosome::recache() {
    update_retro_ = true;
    update_ltr_ = true;
    update_len_ = true;
}

fluke::Chromosome::ce_iter
fluke::Chromosome::upstreamSelect( ce_iter i ) {
    // pre: IsTrueDstream( i )
    ce_riter ii( i );
    ce_riter jj = chro_->rend();
    // BLS for first non-binding site
    while( ii != jj ) {
        if( IsBindingSite()( *ii ) ) {
            ++ii;
        } else {
            jj = ii;
        }
    }
    return jj.base();
}

boost::tuple< fluke::Chromosome*, fluke::Chromosome::ce_iter >
fluke::Chromosome::randChromosomeElement() {
    // pre: chro_ is initialised
    std::list< ce_iter > aux;
    ce_iter i = chro_->begin();
    if( IsTrueDstream()( *i ) || IsBindingSite()( *i ) || IsRepeat()( *i ) ) {
        aux.push_back( i );
    }
    ++i;
    while( i != chro_->end() ) {
        if( IsTrueDstream()( *i ) || IsBindingSite()( *i ) ) {
            aux.push_back( i );
        } else if( IsRepeat()( *i ) && 
                   !IsRetroposon()( *( boost::prior( i ) ) ) ) {
            aux.push_back( i );
        }
        ++i;
    }
    return boost::make_tuple( this, 
        *( random_element( aux.begin(), aux.end(), rand_range< int > ) ) );
}

fluke::Chromosome::tag_container
fluke::Chromosome::essentialTags() const {
    tag_container result;
    if( chro_->empty() ) return result;
    
    for( const_ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        OrdinaryDownstream *aux = dynamic_cast< OrdinaryDownstream* >( *i );
        if( aux /* && IsBindingSite()( *( boost::prior( i ) ) ) */ ) {
            result.push_back( aux->tag() );
        }
    }
    return result;    
}

fluke::Chromosome::tag_container
fluke::Chromosome::moduleTags( int module ) const {
    // return a vector with all the tags in it (incl. duplicates) that have an
    // upstream region
    tag_container result;
    if( chro_->empty() ) return result;
    
    for( const_ce_iter i = chro_->begin(); i != chro_->end(); ++i ) {
        ModuleDownstream *aux = dynamic_cast< ModuleDownstream* >( *i );
        if( aux /* && IsBindingSite()( *( boost::prior( i ) ) ) */ ) {
            if( aux->module() == module ) {
                result.push_back( aux->tag() );
            }
        }
    }
    return result;
}

int
fluke::Chromosome::nrRepeats() const {
    if( update_ltr_ ) {
        update_ltr_ = false;
        nr_ltr_ = std::count_if( chro_->begin(), chro_->end(), IsRepeat() ); 
    }
    return nr_ltr_;
}

int
fluke::Chromosome::nrRetroposons() const {
    if( update_retro_ ) {
        update_retro_ = false;
        nr_retroposons_ = nrRetroposons( chro_->begin(), chro_->end() );
    }
    return nr_retroposons_;
}

int
fluke::Chromosome::nrRetroposons( ce_iter first, ce_iter last ) const {
    return std::count_if( first, last, IsRetroposon() );
}

// Reviewer 2
bool
fluke::Chromosome::oneCentromere() const {
    return std::count_if( chro_->begin(), chro_->end(), IsCentromere() ) == 1;
}

int
fluke::Chromosome::size() const { 
    if( update_len_ ) {
        update_len_ = false;
        len_ = chro_->size();
    }
    return len_;
}

void
fluke::Chromosome::write( std::ostream &os ) const {
    // pre: chro_ is initialised
    os << "<chromosome len=\"" << chro_->size() << "\">\n";
    ce_iter i = chro_->begin(); 
    while( i != chro_->end() ) {
        os << **i;
        ++i;
    }
    // and now the mutation rates
    os << "<rates>\n";
    os << "<copy gene=\"" << cp_gene_rate_ 
       << "\" retro=\"" << cp_tp_rate_ << "\"/>\n";
    os << "<remove gene=\"" << rm_gene_rate_ << "\" retro=\""
       << rm_tp_rate_ << "\" repeat=\"" << rm_ltr_rate_ << "\"/>\n";
    os << "<break dsb=\"" << dsb_recombination_ << "\"/>\n";
    os << "</rates>\n";
    os << "</chromosome>\n";
}

void
fluke::Chromosome::copyRates( const Chromosome &s, Chromosome &t ) const {
    // retroposon
    t.cp_tp_rate_ = s.cp_tp_rate_;
    t.rm_tp_rate_ = s.rm_tp_rate_;
    t.rm_ltr_rate_ = s.rm_ltr_rate_;
    t.new_tp_rate_ = s.new_tp_rate_;
    // bsites
    t.nw_bs_rate_ = s.nw_bs_rate_;
    t.cp_bs_rate_ = s.cp_bs_rate_;
    t.rm_bs_rate_ = s.rm_bs_rate_;
    // gene insdel
    t.cp_gene_rate_ = s.cp_gene_rate_;
    t.rm_gene_rate_ = s.rm_gene_rate_;
    // recombination
    t.dsb_recombination_ = s.dsb_recombination_;
    // and the rates
    t.mut_rate_ = s.mut_rate_;
    t.mut_step_ = s.mut_step_;
    t.dsb_step_ = s.dsb_step_;
    t.retro_step_ = s.retro_step_;
}

void
fluke::Chromosome::copyRetroposonRate( double f )
{ cp_tp_rate_ = f; }

double
fluke::Chromosome::copyRetroposonRate() const
{ return cp_tp_rate_; }

void
fluke::Chromosome::removeRetroposonRate( double f )
{ rm_tp_rate_ = f; }

double
fluke::Chromosome::removeRetroposonRate() const
{ return rm_tp_rate_; }

void
fluke::Chromosome::removeRepeatRate( double f )
{ rm_ltr_rate_ = f; }

double
fluke::Chromosome::removeRepeatRate() const
{ return rm_ltr_rate_; }

void
fluke::Chromosome::newRetroposonRate( double f )
{ new_tp_rate_ = f; }

double
fluke::Chromosome::newRetroposonRate()
{ return new_tp_rate_; }

void 
fluke::Chromosome::recombinationRate( double f )
{ dsb_recombination_ = f; }

double
fluke::Chromosome::recombinationRate() const
{ return dsb_recombination_; }

void
fluke::Chromosome::newBsiteRate( double f )
{ nw_bs_rate_ = f; }

double
fluke::Chromosome::newBsiteRate() const
{ return nw_bs_rate_; }

void
fluke::Chromosome::copyBsiteRate( double f )
{ cp_bs_rate_ = f; }

double
fluke::Chromosome::copyBsiteRate() const
{ return cp_bs_rate_; }

void
fluke::Chromosome::removeBsiteRate( double f )
{ rm_bs_rate_ = f; }

double
fluke::Chromosome::removeBsiteRate() const
{ return rm_bs_rate_; }

void
fluke::Chromosome::copyGeneRate( double f ) 
{ cp_gene_rate_ = f; }

double 
fluke::Chromosome::copyGeneRate() const
{ return cp_gene_rate_; }

void 
fluke::Chromosome::removeGeneRate( double f )
{ rm_gene_rate_ = f; }

double 
fluke::Chromosome::removeGeneRate() const
{ return rm_gene_rate_; }

void
fluke::Chromosome::mutationRate( double f )
{ mut_rate_ = f; }

double
fluke::Chromosome::mutationRate() const
{ return mut_rate_; }

void
fluke::Chromosome::mutationStep( double f )
{ mut_step_ = f; }

double
fluke::Chromosome::mutationStep() const
{ return mut_step_; }

void
fluke::Chromosome::retroStep( double f )
{ retro_step_ = f; }

double
fluke::Chromosome::retroStep() const
{ return retro_step_; }

void
fluke::Chromosome::dsbStep( double f )
{ dsb_step_ = f; }

double
fluke::Chromosome::dsbStep() const
{ return dsb_step_; }

void
fluke::Chromosome::mutationScheme( MutateRates *m )
{ rate_mutator_ = m; }
