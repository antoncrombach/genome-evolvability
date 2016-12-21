//
// Transcription network, the second level.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_SIMPLE_NETWORK_H_
#define _FLUKE_SIMPLE_NETWORK_H_

#include "defs.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "downstream.hh"
#include "bsite.hh"
#include "transfac.hh"

namespace fluke {

    /// \class SimpleNetwork
    /// \brief Boolean transcription network.
    ///
    /// In the future this class will be refactored to be a child of some 
    /// base class. Yet, at this point we have no idea what such a base class
    /// would look like..
    class SimpleNetwork {
        public:
            /// Typedef of a property of a vertex; associating vertices with
            /// downstreams
            typedef boost::property< boost::vertex_name_t, Downstream* > 
                gene_ptr;
            /// Typedef of the graph storage structure
            typedef boost::adjacency_list< boost::vecS, boost::vecS, 
                boost::directedS, gene_ptr > boolean_net;
            /// Typedef of mapping the network to vertex names
            typedef boost::property_map< boolean_net, 
                boost::vertex_name_t >::type gene_ptr_map;
            /// Typedef of reverse mapping, from downstreams to vertices
            typedef std::map< Downstream*, 
                boolean_net::vertex_descriptor > gene_net_map;

        public:
            /// Constructor
            SimpleNetwork();
            /// Constructor given a genome
            SimpleNetwork( const Genome & );
            /// Copy constructor
            SimpleNetwork( const SimpleNetwork & );
            /// Destructor
            ~SimpleNetwork();
            /// Copy given network into \c this
            void copy( const SimpleNetwork & );

            /// Build a network from the given genome
            void build( const Genome & );
            
            /// Write the network in graphViz format
            void writeViz( std::ostream & );
            /// [\b not \b implemented \b yet]
            void readViz( std::istream & );

        private:
            boolean_net *graph_;
    };

    /// \class DownstreamLabelWriter
    /// \brief Write the label of a downstream region. 
    class DownstreamLabelWriter {
        public:
            /// Constructor given a gene to downstream mapping
            DownstreamLabelWriter( SimpleNetwork::gene_ptr_map &gp  ) 
                : lbl_( gp ) {};

            /// Output a label of the downstream
            void operator()( std::ostream &os, 
              const SimpleNetwork::boolean_net::vertex_descriptor &v ) const {
                std::string aux;
                if( lbl_[ v ] != 0 ) {
                    aux = lbl_[ v ]->writeLabel();
                }
                os << "[label=\"" << aux << "\"]"; 
            }

        private:
            SimpleNetwork::gene_ptr_map lbl_;
    };
}
#endif

