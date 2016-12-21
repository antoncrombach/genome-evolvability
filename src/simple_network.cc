//
// Implementation of a boolean transcription network
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "simple_network.hh"

fluke::SimpleNetwork::SimpleNetwork() {
    graph_ = new boolean_net();
}

fluke::SimpleNetwork::SimpleNetwork( const Genome &g ) {
    graph_ = new boolean_net();
    build( g );
}

fluke::SimpleNetwork::SimpleNetwork( const SimpleNetwork &bn ) {
    copy( bn );
}

fluke::SimpleNetwork::~SimpleNetwork() {
    delete graph_;
}

void 
fluke::SimpleNetwork::copy( const SimpleNetwork &bn ) {
    // deep copy?
    graph_ = bn.graph_;
}

void 
fluke::SimpleNetwork::build( const Genome &g ) {
    // First gather all tf's and put them in the graph already.
    gene_net_map tf_map;
    const std::list< Chromosome* > aux = g.chromosomes();
    for( Genome::const_chromos_iter i = aux.begin(); i != aux.end(); ++i ) {
        const std::list< ChromosomeElement* > bux = ( **i ).elements();
        for( Chromosome::const_ce_iter j = bux.begin(); j != bux.end(); ++j ) {
            TranscriptionFactor *tf = 
                dynamic_cast< TranscriptionFactor * >( *j );
            if( tf ) {
                boolean_net::vertex_descriptor vux = 
                    boost::add_vertex( tf, *graph_ );
                tf_map.insert( std::make_pair( tf, vux ) );
            }
        }
    }

    // Second re-iterate the genome and for each bsite check if anything binds
    // to it... naive implementation, should be possible to make it faster!
    for( Genome::const_chromos_iter i = aux.begin(); i != aux.end(); ++i ) {
        const std::list< ChromosomeElement* > bux = ( **i ).elements();
        bool upstream = false;
        boolean_net::vertex_descriptor vux;
        for( Chromosome::const_ce_riter j = bux.rbegin();
                j != bux.rend(); ++j ) {
            if( typeid( **j ) == typeid( BindingSite ) ) {
                if( upstream ) {
                    BindingSite *bs = dynamic_cast< BindingSite* >( *j );
                    gene_net_map::iterator tf = tf_map.begin();
                    while( tf != tf_map.end() ) {
                        TranscriptionFactor *cux = 
                            dynamic_cast< TranscriptionFactor* >( tf->first );
                        if( cux->binds( *bs ) ) {
                            boost::add_edge( tf->second, vux, *graph_ );
                        }
                        ++tf;
                    }
                }
            } else if( typeid( **j ) == typeid( Retroposon ) ) {
                upstream = false;
            } else if( typeid( **j ) == typeid( Repeat ) ) {
                upstream = false;
            } else if( typeid( **j ) == typeid( TranscriptionFactor ) ) {
                upstream = true;
                vux = tf_map[ dynamic_cast< TranscriptionFactor* >( *j ) ];
            } else {
                upstream = true;
                vux = boost::add_vertex( 0, *graph_ );
            }
        }
    }
}

void 
fluke::SimpleNetwork::writeViz( std::ostream &os ) {
    gene_ptr_map aux = boost::get( boost::vertex_name, *graph_ );
    DownstreamLabelWriter bux( aux );
    boost::write_graphviz( os, *graph_, bux );
}

void 
fluke::SimpleNetwork::readViz( std::istream &is ) {
    //boost::read_graphviz( is, *graph_ );
}

