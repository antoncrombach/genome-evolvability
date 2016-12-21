//
// Implementation of a grid/population.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "population.hh"
#include "agent.hh"
#include "genome.hh"
#include "logger.hh"
#include "duo_agent.hh"

bool fluke::Population::shuffle_ = false;
double fluke::Population::threshold_ = 0.0;
int fluke::Population::nr_agent_types_ = 0;
std::string fluke::Population::placement_ = "random";

fluke::Population::Population() 
    : plane_one_(), plane_two_(), 
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(), agent_view_size_( 3 ) {
    scaling_ = new NoScaling();
    selection_ = new RandSelection();
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

fluke::Population::Population( int x, int y, std::vector< Agent* > &vag, 
        ScalingScheme *sca, SelectionScheme *sel ) 
    : plane_one_( boost::extents[ x ][ y ] ), 
      plane_two_( boost::extents[ x ][ y ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(),  
      scaling_( sca ), selection_( sel ), agent_view_size_( 3 ) {
    // pre: x * y > |vag|
    locations( shuffle_locs_ );
    zero( reading );
    zero( writing );
    insert( vag );
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

fluke::Population::Population( int x, int y, std::vector< Agent* > &vag, 
    const std::vector< Location > &loc, 
    ScalingScheme *sca, SelectionScheme *sel ) 
    : plane_one_( boost::extents[ x ][ y ] ),
      plane_two_( boost::extents[ x ][ y ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(),
      scaling_( sca ), selection_( sel ), agent_view_size_( 3 ) {
    // pre: x * y > |vag|
    locations( shuffle_locs_ );
    zero( reading );
    zero( writing );
    insertAt( vag, loc );
    async_agent_obs_ = 0;
    async_env_change_ = 0;
    async_dsbs_ = 0;
}

fluke::Population::Population( const Population &pop ) 
    : plane_one_( 
        boost::extents[ pop.write_grid_->shape()[ 0 ] ][ pop.write_grid_->shape()[ 1 ] ] ),
      plane_two_( 
        boost::extents[ pop.write_grid_->shape()[ 0 ] ][ pop.write_grid_->shape()[ 1 ] ] ),
      write_grid_( &plane_one_ ), read_grid_( &plane_two_ ),
      write_agents_(), read_agents_(), scaling_( 0 ), 
      selection_( 0 ), agent_view_size_( 3 ) {
    // copy deep, all the way down!!
    zero( reading );
    zero( writing );
    // copy agents in map and grid
    for( const_map_ag_iter i = pop.write_agents_.begin(); 
        i != pop.write_agents_.end(); ++i ) {
        // some cloning
        Agent *aux = i->first->clone();
        insertAt( aux, i->second );
    }
    // copy vector of locations
    for( const_loc_iter i = pop.shuffle_locs_.begin(); 
        i != pop.shuffle_locs_.end(); ++i ) {
        Location aux( *i );
        shuffle_locs_.push_back( aux );
    }
    // last but not least...
    scaling_ = pop.scaling_->clone();
    selection_ = pop.selection_->clone();
    model_ = pop.model_;
    async_dsbs_ = 0;
    async_agent_obs_ = 0;
    async_env_change_ = 0;
}

fluke::Population::~Population() {
    // how to deal with the shadow_agents_??
    agents_map aux, bux;
    std::set_intersection( write_agents_.begin(), write_agents_.end(), 
        read_agents_.begin(), read_agents_.end(), 
        std::inserter( aux, aux.begin() ) );
    std::set_symmetric_difference( write_agents_.begin(), write_agents_.end(), 
        read_agents_.begin(), read_agents_.end(), 
        std::inserter( bux, bux.begin() ) );
    // first explicitly delete the agents
    for( const_map_ag_iter i = aux.begin(); i != aux.end(); ++i ) {
        delete i->first;
    }
    for( const_map_ag_iter i = bux.begin(); i != bux.end(); ++i ) {
        delete i->first;
    }
    // and the rest    
    delete scaling_;
    delete selection_;
    if( async_dsbs_ != 0 ) {
        delete async_dsbs_;
    }
    if( async_agent_obs_ != 0 ) {
        delete async_agent_obs_;
    }
    if( async_env_change_ != 0 ) {
        delete async_env_change_;
    }
}

std::vector< fluke::Agent* >
fluke::Population::moore( const Location &loc ) {
    // torus border
    int n = read_grid_->shape()[ 0 ];
    int m = read_grid_->shape()[ 1 ];
    // do somethin smart...
    std::vector< Agent* > result;
    int i = (loc.x - 1 + n) % n;
    while( i != (loc.x + 2) % n ) {
        int j = (loc.y - 1 + m) % m;
        while( j != (loc.y + 2) % m ) {
            result.push_back( ( *read_grid_ )[ i ][ j ] );
            j = (j + 1) % m;
        }
        i = (i + 1) % n;
    }
    return result;
}

std::vector< fluke::Agent* >
fluke::Population::moore( const Agent &ag ) {
    return moore( read_agents_[ const_cast< Agent* >( &ag ) ] );
}

fluke::Location
fluke::Population::whereIs( const Agent &ag ) {
    return read_agents_[ const_cast< Agent* >( &ag ) ];
}

void
fluke::Population::initialise() {
    // pre: model_ is set
#ifdef DEBUG
    cout << "! init population" << endl;
#endif
    for( map_ag_iter i = write_agents_.begin(); 
        i != write_agents_.end(); ++i ) {
        AgentTag aux = AgentTag( model_->now(), i->second.x, i->second.y, 0 );
        i->first->myTag( aux );
    }
#ifdef DEBUG
    cout << "! population ready" << endl;
    cout << "  # agent types: " << nr_agent_types_ << "\n"
         << "  placement    : " << placement_ << endl;
#endif
    // do a first shuffle and update of the planes
    swap();
    if( shuffle_ ) shuffle();
}

void 
fluke::Population::evaluate( const Environment &env ) {
    // first dump population in file
    // hack!
    if( async_env_change_ != 0 ) {
        async_env_change_->update( this );
    }
    // end hack!
    for( map_ag_iter i = write_agents_.begin(); 
        i != write_agents_.end(); ++i ) {
        i->first->evaluate( env );
    }
}

void 
fluke::Population::step() {
    // deterministic synchronous stepping
    // visit every site
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            if( ( *read_grid_ )[ i ][ j ] != 0 ) {
                // occupied spot
                ( *read_grid_ )[ i ][ j ]->step( *this );
                if( ( *read_grid_ )[ i ][ j ]->dying() ) {
                    shallowErase( ( *write_grid_ )[ i ][ j ] );
                }
            } else {
                // empty spot
                // get nbh
                Location nux;
                nux.x = i;
                nux.y = j;
                std::vector< Agent* > aux( moore( nux ) );
                // remove non-agents
                aux.erase( std::remove_if( aux.begin(), aux.end(), 
                        IsNotAvailable() ), aux.end() );
                if( aux.size() > 0 ) {
                    // scores of agents in nbh
                    std::vector< double > cux;
                    for( ag_iter k = aux.begin(); k < aux.end(); ++k ) {
                        cux.push_back( ( **k ).score() );
                    }
                    // scale the scores
                    scaling_->scale( cux );
                    // check for sum of fitness
                    double bux = threshold_ -
                        std::accumulate( cux.begin(), cux.end(), 0.0 );
                    if( bux > 0.0 ) {
                        aux.push_back( 0 );
                        cux.push_back( bux );
                    }
                    // select an agent (or a null pointer)
                    Agent *eux = selection_->select( aux, cux );
                    if( eux != 0 ) {
                        // log this agent be4 mutations
                        if( async_agent_obs_ != 0 ) {
                            async_agent_obs_->update( eux );
                        }
                        // spawn a sibling
                        Agent *fux = eux->sibling();
                        // get the mother tag and set ancestor tags of children
                        AgentTag gux = eux->myTag();
                        eux->parentTag( gux );
                        fux->parentTag( gux );
                        // update children's tag
                        if( gux.time == model_->now() ) {
                            // increase index, agent is part of a 'cascade'
                            ++gux.i;
                        }
                        gux.time = model_->now();
                        eux->myTag( gux );
                        fux->myTag( AgentTag( model_->now(), i, j, 0 ) );
                        // evaluate both in the environment
                        eux->evaluate( model_->environment() );
                        fux->evaluate( model_->environment() );
                        // and insert it in the grid
                        insertAt( fux, nux );
                        // log after mutations what happened (dsbs)
                        if( async_dsbs_ != 0 ) {
                            DuoAgent hux( eux, fux );
                            async_dsbs_->update( &hux );
                        }
                    }
                }
            }
        }
    }
    
    // keep everything consistent
    swap();
    if( shuffle_ ) shuffle();    
}

void
fluke::Population::finish() {
    // log all agents to the ancestor tracing observer
    if( async_agent_obs_ != 0 ) {
        for( map_ag_iter i = write_agents_.begin(); 
            i != write_agents_.end(); ++i ) {
            async_agent_obs_->update( i->first );
        }
    }
}

void 
fluke::Population::insert( std::vector< Agent* > &ag ) {
    // insert at locations according to placement
    if( placement_ == "patch" ) {
        // assuming 2 different agent types, just sorting them
        std::sort( ag.begin(), ag.end(), SmallerTypeThan() );
        // and need to know how much of each
        ag_iter tt = std::lower_bound( ag.begin(), ag.end(), 
            ag.back(), SmallerTypeThan() );
        std::vector< int > aux;
        aux.push_back( std::distance( ag.begin(), tt ) );
        aux.push_back( std::distance( tt, ag.end() ) );
        // get the patches
        std::vector< std::vector< Location > > bux = patchEmptyLocations( aux );
        // only supporting two patches
        insertAt( ag.begin(), tt, bux[ 0 ] );
        insertAt( tt, ag.end(), bux[ 1 ] );
    } else {
        // spread over entire field
        std::vector< Location > aux = randEmptyLocations( ag.size() );
        insertAt( ag, aux );
    }
}

void
fluke::Population::insertAt( ag_iter first, ag_iter last, 
        const std::vector< Location > &loc ) {
    ag_iter i = first;
    const_loc_iter j = loc.begin();
    while( i != last ) {
        insertAt( *(i++), *(j++) );
    }
}    

void 
fluke::Population::insertAt( 
        std::vector< Agent* > &ag, const std::vector< Location > &loc ) {
    // pre: ag.size() == loc.size()
    // pos: memory management of agents is responsibility of this class
    ag_iter i = ag.begin();
    const_loc_iter j = loc.begin();
    while( i != ag.end() ) {
        insertAt( *(i++), *(j++) );
    }
}

void
fluke::Population::insertAt( Agent* ag, const Location &loc ) {
    write_agents_[ ag ] = loc;
    ( *write_grid_ )[ loc.x ][ loc.y ] = ag;
}

void 
fluke::Population::erase( std::vector< Agent* > &ag ) {
    std::vector< Location > result;
    for( ag_iter i = ag.begin(); i != ag.end(); ++i ) {
        result.push_back( write_agents_[ *i ] );
    }
    eraseAt( result );
}

void
fluke::Population::erase( Agent* ag ) {
    eraseAt( write_agents_[ ag ] );
}

void 
fluke::Population::eraseAt( const std::vector< Location > &loc ) {
    for( const_loc_iter i = loc.begin(); i != loc.end(); ++i ) {
        eraseAt( *i );
    }
}

void
fluke::Population::eraseAt( Location loc ) {
    write_agents_.erase( ( *write_grid_ )[ loc.x ][ loc.y ] );
    delete ( *write_grid_ )[ loc.x ][ loc.y ];
    ( *write_grid_ )[ loc.x ][ loc.y ] = 0;
}

void
fluke::Population::shallowErase( Agent* ag ) {
    shallowEraseAt( write_agents_[ ag ] );
}

void
fluke::Population::shallowEraseAt( Location loc ) {
    write_agents_.erase( ( *write_grid_ )[ loc.x ][ loc.y ] );
    ( *write_grid_ )[ loc.x ][ loc.y ] = 0;
}

std::vector< fluke::Location >
fluke::Population::randEmptyLocations( int n ) {
    // find at most n random empty locations
    // note: if necessary, it can be optimised. Yet we're not planning on
    // putting our energy to the optimisation of this method. It's only
    // used once during initialisation of the program FIX
    std::vector< Location > result;
    locations( result );
    std::random_shuffle( result.begin(), result.end(), rand_range< int > );
    loc_iter i( result.begin() );
    loc_iter j( boost::prior( result.end() ) );
    while( i != j ) {
        if( ( *write_grid_ )[ i->x ][ i->y ] == 0 ) {
            ++i;
        } else {
            std::iter_swap( i, j );
            --j;
        }
    }

    int k = std::distance( result.begin(), i );
    if( k < n ) n = k;
    result.erase( result.begin() + n, result.end() );
    return result;
}

std::vector< std::vector< fluke::Location > >
fluke::Population::patchEmptyLocations( std::vector< int > n ) {
    // find n random locations near each other
    // no general algorithm... only nr_agent_types = {1, 2} supported
    std::vector< std::vector< Location > > result;
    std::vector< Location > aux;
    locations( aux );
    std::random_shuffle( aux.begin(), aux.end(), rand_range< int > );
    int bux = static_cast< int >( write_grid_->shape()[ 0 ] ) / 2;
    // do it
    for( std::vector< int >::iterator ii = n.begin(); ii != n.end(); ++ii ) {
        loc_iter jj = aux.begin();
        result.push_back( std::vector< Location >() );
        int k = 0; 
        while( k != *ii ) {
            if( jj->x >= std::distance( n.begin(), ii ) * bux &&
                jj->x < ( std::distance( n.begin(), ii ) + 1 ) * bux ) {
                result.back().push_back( *jj );
                ++k;
            }
            // no check if jj goes past end!!
            ++jj;
        }
    }
    return result;
}

void
fluke::Population::locations( std::vector< Location > &ll ) {
    // grids are of same size
    ll.clear();
    ll.reserve( read_grid_->shape()[ 0 ] * read_grid_->shape()[ 1 ] );
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            Location aux;
            aux.x = i;
            aux.y = j;
            ll.push_back( aux );
        }
    }
}

void
fluke::Population::zero( grid_type tt ) {
    agents_grid *g = write_grid_;
    if( tt == reading ) {
        g = read_grid_;
    }
    
    for( uint i = 0; i < g->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < g->shape()[ 1 ]; ++j ) {
            ( *g )[ i ][ j ] = static_cast< Agent* >( 0 );
        }
    }
}

void
fluke::Population::swap() {
    // swap grids and mapping, shadow is readable now...
    std::swap( write_grid_, read_grid_ );
    write_agents_.swap( read_agents_ );
    
    // keep data consistent
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            if( (* read_grid_ )[ i ][ j ] != 0 && 
                ( *write_grid_ )[ i ][ j ] == 0 ) {
                // insertion happened last time
                Location nux( i, j );
                insertAt( ( *read_grid_ )[ i ][ j ], nux ); 
            } else if( (* read_grid_ )[ i ][ j ] == 0 &&
                ( *write_grid_ )[ i ][ j ] != 0 ) {
                // deletion happened last time
                Location nux( i, j );
                eraseAt( nux );
            }
        }
    }
}

void
fluke::Population::swap( grid_type tt, const Location &a, const Location &b ) {
    agents_grid *g = write_grid_;
    agents_map *m = &write_agents_;
    if( tt == reading ) {
        g = read_grid_;
        m = &read_agents_;
    }
    // first swap in map
    if( ( *g )[ a.x ][ a.y ] != 0 ) {
        ( *m )[ ( *g )[ a.x ][ a.y ] ] = b;
    }
    if( ( *g )[ b.x ][ b.y ] != 0 ) {
        ( *m )[ ( *g )[ b.x ][ b.y ] ] = a;
    }
    // now swap in grid
    std::swap( ( *g )[ a.x ][ a.y ], ( *g )[ b.x ][ b.y ] );
}

void
fluke::Population::shuffle() {
    // assuming grids are consistent
    std::random_shuffle( shuffle_locs_.begin(), shuffle_locs_.end(), 
        rand_range< int > );
#ifdef DEBUG
    for( loc_iter i = shuffle_locs_.begin(); i != shuffle_locs_.end(); ++i ) {
        cout << i->x << ", " << i->y << ": ";
    }
    cout << endl;
#endif
    const_loc_iter aux = shuffle_locs_.begin();
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            Location bux( i, j );
            swap( reading, *aux, bux );
            swap( writing, *aux, bux );
            ++aux;
        }
    }
}

void
fluke::Population::attach1( AsyncLogObserver *l ) {
    async_agent_obs_ = l;
}

void
fluke::Population::attach2( AsyncLogObserver *l ) {
    async_env_change_ = l;
}

void
fluke::Population::attach3( AsyncLogObserver *l ) {
    async_dsbs_ = l;
}

void
fluke::Population::closeAll() {
    if( async_agent_obs_ != 0 ) {
        async_agent_obs_->closeLog();
    }
    if( async_env_change_ != 0 ) {
        async_env_change_->closeLog();
    }
    if( async_dsbs_ != 0 ) {
        async_dsbs_->closeLog();
    }
}

void
fluke::Population::threshold( double f ) 
{ threshold_ = f; }
    
double
fluke::Population::threshold()
{ return threshold_; }

void
fluke::Population::nrAgentTypes( int t )
{ nr_agent_types_ = t; }

int
fluke::Population::nrAgentTypes()
{ return nr_agent_types_; }

void
fluke::Population::placement( const std::string &s )
{ placement_ = s; }

std::string
fluke::Population::placement()
{ return placement_; }

void
fluke::Population::shuffling( bool t )
{ shuffle_ = t; }

bool
fluke::Population::shuffling()
{ return shuffle_; }

bool
fluke::Population::hasEveryAgentType() const {
    std::vector< uint > aux( nr_agent_types_, 0 );
    const_map_ag_iter i = write_agents_.begin();
    const_map_ag_iter j = write_agents_.end();
/*#ifdef DEBUG
    std::cout << "! has all agents\n";
#endif*/
    while( i != j ) {
/*#ifdef DEBUG
    std::cout << i->first->type() << ": " 
              << i->second.x << ", " << i->second.y << std::endl;
#endif*/
        if( std::accumulate( aux.begin(), aux.end(), 0 ) != nr_agent_types_ ) {
            aux[ i->first->type() - 1 ] = 1;
            ++i;
        } else {
            j = i;
        }
    }
    return std::accumulate( aux.begin(), aux.end(), 0 ) == nr_agent_types_;
}

void
fluke::Population::write( std::ostream &os ) const {
    // write the 'shadow' plane
    os << "Shadow/read plane:\n";
    for( uint i = 0; i < read_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < read_grid_->shape()[ 1 ]; ++j ) {
            os.width( 10 );
            if( (* read_grid_ )[ i ][ j ] != 0 ) {
                os << (* read_grid_ )[ i ][ j ];
            } else {
                os << ".";
            }
            os << ' ';
        }
        os << "\n";
    }
    os << "Shadow/read map:\n";
    for( const_map_ag_iter i = read_agents_.begin();
        i != read_agents_.end(); ++i ) {
            os << i->first << "\n";
    }
    // write the 'readable' plane
    os << "Current/write plane:\n";
    for( uint i = 0; i < write_grid_->shape()[ 0 ]; ++i ) {
        for( uint j = 0; j < write_grid_->shape()[ 1 ]; ++j ) {
            os.width( 10 );
            if( (* write_grid_ )[ i ][ j ] != 0 ) {
                os << (* write_grid_ )[ i ][ j ];
            } else {
                os << ".";
            }
            os << ' ';
        }
        os << "\n";
    }
    os << "Current/write map:\n";
    for( const_map_ag_iter i = write_agents_.begin();
        i != write_agents_.end(); ++i ) {
            os << i->first << "\n";
    }

}
