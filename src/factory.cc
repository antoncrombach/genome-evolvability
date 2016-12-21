//
// Implementation of the model factory
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "factory.hh"
#include "bsite.hh"
#include "transfac.hh"
#include "module_dstream.hh"
#include "ordinary_dstream.hh"
#include "retroposon.hh"
#include "repeat.hh"
#include "chromosome.hh"
#include "genome.hh"
#include "scaling.hh"
#include "selection.hh"
#include "simple_agent.hh"
#include "module_agent.hh"
#include "population.hh"
#include "environment.hh"
#include "agent_reader.hh"
#include "population_reader.hh"
#include "mutate_rates.hh"
// reviewer 2
#include "centromere.hh"

XERCES_CPP_NAMESPACE_USE

fluke::Factory::Factory() {
    // empty managers
    bs_mng_ = 0;
    tf_mng_ = 0;
    downstream_tag_ = 0;
}

fluke::Factory::Factory( Config *gc ) : conf_( gc ) {
    // empty managers
    bs_mng_ = 0;
    tf_mng_ = 0;
    downstream_tag_ = 0;

}

fluke::Factory::~Factory() {
    if( bs_mng_ != 0 )
        delete bs_mng_;
    if( tf_mng_ != 0 )
        delete tf_mng_;
}

fluke::BindingSite* 
fluke::Factory::bsite() {
    return new BindingSite( bsiteManager()->allocRandShortSeq() );
}

fluke::TranscriptionFactor* 
fluke::Factory::transfac() {
    return new TranscriptionFactor( tag(), 
            transfacManager()->allocRandShortSeq() );
}

fluke::ModuleDownstream*
fluke::Factory::moduleDstream( int m ) {
    return new ModuleDownstream( tag(), m );
}

fluke::OrdinaryDownstream*
fluke::Factory::ordinaryDstream() {
    return new OrdinaryDownstream( tag() );
}

fluke::Retroposon* 
fluke::Factory::retroposon() {
    return new Retroposon( tag() );
}

fluke::Repeat*
fluke::Factory::repeat() {
    return new Repeat();
}

fluke::Centromere*
fluke::Factory::centromere() {
    return new Centromere();
}

fluke::Chromosome* 
fluke::Factory::chromosome( int type ) {
    // auxilary typedef
    std::list< ChromosomeElement* > *ll = 
        organiseChromosome( chromosomeParts( type ), 
            conf_->optionAsDouble( "organised", type ) );
    
    // scan list and flank all retroposons with repeats
    // scan list and add bsites to the downstreams
    int nr_bs = conf_->optionAsInt( "bsites", type );
    Chromosome::ce_iter i = ll->begin();
    while( i != ll->end() ) {
        if( Chromosome::IsRetroposon()( *i ) ) {
            ll->insert( i, repeat() );
            ll->insert( boost::next( i ), repeat() );
        } else if( Chromosome::IsTrueDstream()( *i ) ) {
            for( int j = 0; j != nr_bs; ++j ) {
                ll->insert( i, bsite() );
            }
        }
        ++i;
    }
    
    // init with right values for this type of agent
    Chromosome *result = new Chromosome( 0, ll );
    result->newBsiteRate( conf_->optionAsDouble( "new_bsite", type ) );
    result->copyBsiteRate( conf_->optionAsDouble( "cp_bsite", type ) );
    result->removeBsiteRate( conf_->optionAsDouble( "rm_bsite", type ) );
    result->copyGeneRate( conf_->optionAsDouble( "cp_gene", type ) );
    result->removeGeneRate( conf_->optionAsDouble( "rm_gene", type ) );
    result->copyRetroposonRate( conf_->optionAsDouble( "cp_tp", type ) );
    result->removeRetroposonRate( conf_->optionAsDouble( "rm_tp", type ) );
    result->newRetroposonRate( conf_->optionAsDouble( "new_tp", type ) );
    result->recombinationRate( 
        conf_->optionAsDouble( "dsb_recombination", type ) );
    result->removeRepeatRate( conf_->optionAsDouble( "rm_ltr", type ) );
    result->mutationStep( conf_->optionAsDouble( "mut_step", type ) );
    result->dsbStep( conf_->optionAsDouble( "dsb_step", type ) );
    result->retroStep( conf_->optionAsDouble( "retro_step", type ) );
    result->mutationRate( conf_->optionAsDouble( "mut_rate", type ) );
    result->mutationScheme( mutateRates( type ) );
    return result;
}

fluke::Genome* 
fluke::Factory::genome( int type ) {
    resetTag();
    // #chromosomes
    int n_ch = conf_->optionAsInt( "chromos", type );

    std::list< Chromosome* > *ll = new std::list< Chromosome* >();
    for( int i = 0; i < n_ch; ++i ) {
        ll->push_back( this->chromosome( type ) );
    }
    return new Genome( ll );   
}

fluke::Agent*
fluke::Factory::agent( int type ) {
    // build agent according to type    
    Agent *result;
    if( conf_->optionAsString( "agent", type ) == "simple" ) {
        result = new SimpleAgent();
    } else if( conf_->optionAsString( "agent", type ) == "module" ) {
        result = new ModuleAgent( type, genome( type ) );
    } else {
        // try to open it as a file
        result = readAgent( conf_->optionAsString( "agent", type ), type );
    }
    result->initialise();
    return result;
}

fluke::Population*
fluke::Factory::population() {
    Population *result = 0;
    // read some population level parameters
    Population::nrAgentTypes( conf_->optionAsInt( "nr_agent_type" ) );
    Population::placement( conf_->optionAsString( "agent_placement" ) );
    Population::shuffling( conf_->optionAsString( "shuffle" ) == "true" );
    Population::threshold( conf_->optionAsDouble( "sum_fitness_threshold" ) );
    // and per agent type stuff
    readAgentConfigurations();
    
    uint xx = conf_->optionAsInt( "grid_x" );
    uint yy = conf_->optionAsInt( "grid_y" );
    bool first_pop = conf_->hasOption( "population_one" );
    bool second_pop = conf_->hasOption( "population_two" );
    bool rand_pop = conf_->optionAsString( "population_start" ) == "random"; 
    if( !first_pop && !second_pop ) {
        std::cout << "Creating one population" << std::endl;
        std::vector< Agent* > aux;
        int cux = conf_->optionAsInt( "nr_agent_type" ) + 1;
        // random start or homogeneous start
        if( rand_pop ) {
            for( int i = 1; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                for( int j = 0; j != bux; ++j ) {
                    aux.push_back( agent( i ) );
                }
            }
        } else {
            for( int i = 1; i != cux; ++i ) {
                int bux = conf_->optionAsInt( "init_nr_agents", i );
                aux.push_back( agent( i ) );
                for( int j = 1; j < bux; ++j ) {
                    aux.push_back( aux.front()->clone() );
                }
            }
        }
        result = new Population( xx, yy, aux,
            scalingScheme(), selectionScheme() );
    } else if( first_pop && !second_pop ) {
        std::cout << "Reading in one population (1)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< Location > loc;
        boost::tie( aux, loc ) = 
            readPopulation( conf_->optionAsString( "population_one" ), 1 );
        // set agent type to 1 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 1 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else if( !first_pop && second_pop ) {
        std::cout << "Reading in one population (2)" << std::endl;
        std::vector< Agent* > aux;
        std::vector< Location > loc;
        boost::tie( aux, loc ) = 
            readPopulation( conf_->optionAsString( "population_one" ), 2 );
        // set agent type to 1 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            //( **ii ).type( 2 );
            ( **ii ).initialise();
        } 
        result = new Population( xx, yy, aux, loc,
            scalingScheme(), selectionScheme() );
    } else { // both populations
        std::cout << "Reading in two populations" << std::endl;
        std::vector< Agent* > aux;
        std::vector< Location > loc;
        boost::tie( aux, loc ) = 
            readPopulation( conf_->optionAsString( "population_one" ), 1 );
        // set agent type to 1 and init
        for( std::vector< Agent* >::iterator ii = aux.begin(); 
            ii != aux.end(); ++ii ) {
            ( **ii ).type( 1 );
            ( **ii ).initialise();
        } 
        // 2nd population
        std::vector< Agent* > bux;
        std::vector< Location > bloc;
        boost::tie( bux, bloc ) = 
            readPopulation( conf_->optionAsString( "population_two" ), 2 );
        // translation of locations
        uint dx = conf_->optionAsInt( "grid_x" ) / 2;
        for( std::vector< Location >::iterator ii = bloc.begin();
            ii != bloc.end(); ++ii ) {
            ii->x += dx;
        }
        // set agent type to 2 and init
        for( std::vector< Agent* >::iterator ii = bux.begin(); 
            ii != bux.end(); ++ii ) {
            ( **ii ).type( 2 );
            ( **ii ).initialise();
        }
        result = new Population( xx, yy, bux, bloc,
            scalingScheme(), selectionScheme() );
    }
    return result;
}

fluke::Environment*
fluke::Factory::environment() {
    // FIX using magic numbers for the number of states
    Environment *result;
    if( conf_->optionAsString( "environment" ) == "constant" ) {
        result = new ConstantEnvironment( 2 );
        ConstantEnvironment *aux =
            dynamic_cast< ConstantEnvironment * >( result );
        // Set modules
        aux->expectedCopies( 0, conf_->optionAsInt( "low_module_a" ) );
        aux->expectedCopies( 1, conf_->optionAsInt( "low_module_b" ) );
    } else if( conf_->optionAsString( "environment" ) == "poisson" ) {
        result = new PoissonEnvironment( 2 );
        PoissonEnvironment *aux =
            dynamic_cast< PoissonEnvironment * >( result );
        // Set one module
        aux->states( 0, conf_->optionAsInt( "low_module_a" ),
                              conf_->optionAsInt( "high_module_a" ) );
        aux->lambda( 0, conf_->optionAsDouble( "lambda_module_a" ) );
        // Set another module
        aux->states( 1, conf_->optionAsInt( "low_module_b" ),
                              conf_->optionAsInt( "high_module_b" ) );
        aux->lambda( 1, conf_->optionAsDouble( "lambda_module_b" ) );       
    } else if( conf_->optionAsString( "environment" ) == "periodic" ) {
        result = new PeriodicEnvironment( 2 );
        PeriodicEnvironment *aux =
            dynamic_cast< PeriodicEnvironment * >( result );
        // Set one module
        aux->states( 0, conf_->optionAsInt( "low_module_a" ),
                              conf_->optionAsInt( "high_module_a" ) );
        aux->period( 0, static_cast< long >( 
                conf_->optionAsDouble( "lambda_module_a" ) ) );
        aux->offset( 0, conf_->optionAsInt( "offset_a" ) );
        // Set another module
        aux->states( 1, conf_->optionAsInt( "low_module_b" ),
                              conf_->optionAsInt( "high_module_b" ) );
        aux->period( 1, static_cast< long >( 
                conf_->optionAsDouble( "lambda_module_b" ) ) );
        aux->offset( 1, conf_->optionAsInt( "offset_b" ) );
    } else {
        result = new ConstantEnvironment( 2 );
    }
    return result;
}

fluke::MutateRates*
fluke::Factory::mutateRates( int type ) {
    MutateRates *result;
    if( conf_->optionAsString( "mutate_scheme", type ) == "fixed" ) {
        result = new FixedMutateRates();
    } else if( conf_->optionAsString( "mutate_scheme", type ) == "linear" ) {
        result = new LinearMutateRates();
    } else if( conf_->optionAsString( "mutate_scheme", type ) == "uniform" ) {
        result = new UniformMutateRates( 
            conf_->optionAsDouble( "uniform_low", type ),
            conf_->optionAsDouble( "uniform_high", type ) );
    } else {
        // default
        result = new FixedMutateRates();
    }
    return result;
} 

fluke::ScalingScheme*
fluke::Factory::scalingScheme() {
    ScalingScheme *result;
    if( conf_->optionAsString( "scaling_scheme" ) == "none" ) {
        result = new NoScaling();
    } else if( conf_->optionAsString( "scaling_scheme" ) == "linear" ) {
        result = new LinearScaling( conf_->optionAsDouble( "base_score" ) );
    } else if( conf_->optionAsString( "scaling_scheme" ) == "power" ) {
        result = new PowerScaling( conf_->optionAsDouble( "base_score" ), 10.0 );
    } else {
        // default
        result = new NoScaling();
    }
    return result;
}

fluke::SelectionScheme*
fluke::Factory::selectionScheme() {
    SelectionScheme *result;
    if( conf_->optionAsString( "selection_scheme" ) == "random" ) {
        result = new RandSelection();
    } else if( conf_->optionAsString( "selection_scheme" ) == "probalistic" ) {
        result = new ProbalisticSelection();
    } else {
        // default
        result = new RandSelection();
    }
    return result;
}

fluke::ShortSeqManager* 
fluke::Factory::bsiteManager() {
    if( bs_mng_ == 0 ) {
        bs_mng_ = new ShortSeqManager( conf_->optionAsString( "alphabet" ),
            conf_->optionAsInt( "length" ), 
            conf_->optionAsDouble( "point_mut_bsite" ) );
        BindingSite::shortSeqManager( bs_mng_ );
    }
    return bs_mng_;
}

fluke::ShortSeqManager* 
fluke::Factory::transfacManager() {
    if( tf_mng_ == 0 ) {
        tf_mng_ = new ShortSeqManager( conf_->optionAsString( "alphabet" ),
            conf_->optionAsInt( "length" ), 
            conf_->optionAsDouble( "point_mut_dstream" ) );
        TranscriptionFactor::shortSeqManager( tf_mng_ );
    }
    return tf_mng_;
}

fluke::Agent*
fluke::Factory::readAgent( std::string fname, int tt ) {
    Agent *result = 0;

    try {
        XMLPlatformUtils::Initialize();
    }
    catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n"
             << message << "\n";
        XMLString::release( &message );
    }

    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );

    AgentReader *doc_handler = new AgentReader( this );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    // set agent type
    doc_handler->type( tt );    
    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Badass lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent
    if( !doc_handler->sawErrors() ) {
        result = doc_handler->agent();
    }
    doc_handler->resetDocument();

    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();
    
    return result;
}

boost::tuple< std::vector< fluke::Agent* >, std::vector< fluke::Location > >
fluke::Factory::readPopulation( std::string fname, int tt ) {
    std::vector< Agent* > result;
    std::vector< Location > loc;

    try {
        XMLPlatformUtils::Initialize();
    }
    catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Error during initialization! :\n"
             << message << "\n";
        XMLString::release( &message );
    }
            
    SAX2XMLReader *parser = XMLReaderFactory::createXMLReader();
    parser->setFeature( XMLUni::fgSAX2CoreValidation, false );
    parser->setFeature( XMLUni::fgSAX2CoreNameSpaces, false );

    PopulationReader *doc_handler = new PopulationReader( this );
    ErrorHandler *error_handler = (ErrorHandler*) doc_handler;
    parser->setContentHandler( doc_handler );
    parser->setErrorHandler( error_handler );    

    // set agent type..
    // not giving type, coz it's read from file
    //doc_handler->type( tt );
    try {
        parser->parse( fname.c_str() );
    } catch( const XMLException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release( &message );
    } catch( const SAXParseException& toCatch ) {
        char* message = XMLString::transcode( toCatch.getMessage() );
        std::cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release( &message );
    } catch( boost::bad_lexical_cast &toCatch ) {
        std::cout << "Badass lexical cast: " << toCatch.what() << "\n";
    } catch( ... ) {
        std::cout << "Unexpected Exception \n" ;
    }

    // get that agent
    if( !doc_handler->sawErrors() ) {
        boost::tie( result, loc ) = doc_handler->population();
    }
    doc_handler->resetDocument();

    // delete parser before calling Terminate
    delete parser;
    delete doc_handler;
    // and call Terminate
    XMLPlatformUtils::Terminate();
    
    return boost::make_tuple( result, loc );
}

void
fluke::Factory::readAgentConfigurations() {
    for( int i = 1; i != conf_->optionAsInt( "nr_agent_type" ) + 1; ++i ) {
        std::string fname =
            "agent" + boost::lexical_cast< std::string >( i ) + ".cfg";
        conf_->parseAgentFile( fname );
        setConfiguration( i );
    }
}

void
fluke::Factory::setConfiguration( int type ) {
    // pre: conf_ is initialised
    ShortSeqManager::maxDistance( conf_->optionAsInt( "max_hamming", type ) );

    SimpleAgent::birthRate( conf_->optionAsDouble( "birth_rate", type ) );
    SimpleAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );

    ModuleAgent::birthRate( conf_->optionAsDouble( "birth_rate", type ) );
    ModuleAgent::deathRate( conf_->optionAsDouble( "death_rate", type ) );
    ModuleAgent::maxDistance( conf_->optionAsInt( "max_distance", type ) );
    ModuleAgent::maxGenomeSize( conf_->optionAsInt( "max_genome_size", type ) );
    ModuleAgent::maxTposons( conf_->optionAsInt( "max_tposons", type ) );
    ModuleAgent::genomePenaltyRate( 
        conf_->optionAsDouble( "genome_size_penalty", type ) );
    ModuleAgent::retroposonPenaltyRate(
        conf_->optionAsDouble( "tposons_penalty", type ) );
}

std::vector< std::list< fluke::ChromosomeElement* > > *
fluke::Factory::chromosomeParts( int type ) {
    // # (essential) dstreams
    int n_ds = conf_->optionAsInt( "dstreams", type );
    // # nr of modules
    int n_mod = conf_->optionAsInt( "modules", type );
    // # nr of genes per module
    int n_ds_mod = conf_->optionAsInt( "dstreams_mod", type );
    // # retroposons ( inbetween other downstreams )
    int n_tp = conf_->optionAsInt( "tposons", type );
    // # single repeat elements ( inbetween downstreams etc )
    int n_rp = conf_->optionAsInt( "repeats", type );
    
    // create chromosome parts
    std::vector< std::list< ChromosomeElement* > > *result = 
        new std::vector< std::list< ChromosomeElement* > >();
    // create (essential) downstreams
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_ds; ++i ) {
        result->back().push_back( ordinaryDstream() );
    }
    // create module downstreams
    for( int i = 0; i < n_mod; ++i ) {
        result->push_back( std::list< ChromosomeElement* >() );
        for( int j = n_ds + i * n_ds_mod; j < n_ds + (i+1) * n_ds_mod; ++j ) {
            result->back().push_back( moduleDstream( i ) );
        }
    }
    // create retroposons
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_tp; ++i ) {
        result->back().push_back( retroposon() );
    }
    // create single repeats
    result->push_back( std::list< ChromosomeElement* >() );
    for( int i = 0; i < n_rp; ++i ) {
        result->back().push_back( repeat() );
    }
    // and add a single centromere (to single repeat elements...)
    result->back().push_back( centromere() );
    return result;
}

std::list< fluke::ChromosomeElement* > *
fluke::Factory::organiseChromosome( 
    std::vector< std::list< ChromosomeElement* > > *parts, double organise ) {
    typedef std::vector< std::list< ChromosomeElement* > >::iterator aux_iter;
    std::list< ChromosomeElement* > *result =
        new std::list< ChromosomeElement* >();
    
    // paste perfect & shuffle a bit
    // last 2 elements of 'parts' are retroposon and repeat elements
    std::list< ChromosomeElement* > ltr = parts->back();
    parts->pop_back();
    std::list< ChromosomeElement* > rp = parts->back();
    parts->pop_back();
    // add repeats
    int aux = 2 * parts->size();
    Chromosome::ce_iter i = ltr.begin();
    while( i != ltr.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // and retroposons
    i = rp.begin();
    while( i != rp.end() ) {
        int bux = static_cast< int >( aux * uniform() );
        if( bux % 2 == 0 ) {
            ( *parts )[ bux / 2 ].push_back( *i );
        } else {
            ( *parts )[ bux / 2 ].push_front( *i );
        }
        ++i;
    }
    // glue together
    for( aux_iter i = parts->begin(); i != parts->end(); ++i ) {
        result->splice( result->end(), *i );
    }
    // shuffle a little bit...
    aux = result->size();
    int nr_elem = static_cast< int >( 0.5 + ( 1.0 - organise ) * aux );
    if(  nr_elem == aux ) {
        // most annoying, copy to vector, from vector
        std::vector< ChromosomeElement* > bux( result->begin(), result->end() );
        std::random_shuffle( bux.begin(), bux.end(), rand_range< int > );
        std::copy( bux.begin(), bux.end(), result->begin() );
    } else {
        // naive implementation
        // take out elements
        std::list< ChromosomeElement* > bux;
        for( int i = 0; i < nr_elem; ++i ) {
            Chromosome::ce_iter cux = boost::next( result->begin(), 
                static_cast< int >( aux * uniform() ) );
            bux.push_back( *cux );
            result->erase( cux );
            --aux;
        }
        // put them back 
        while( !bux.empty() ) {
            result->insert( 
                boost::next( result->begin(), 
                    static_cast< int >( aux * uniform() ) ), bux.back() );
            bux.pop_back();
            ++aux;
        }
    }
    delete parts;
    return result;
}
