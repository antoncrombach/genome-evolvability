//
// The heart of the program 'fluke'. A very simple heart :)
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "fluke.hh"
#include "config.hh"
#include "model.hh"
#include "stream_manager.hh"


fluke::Fluke::Fluke( int argc, char **argv ) 
    : stream_( 0 ), model_( 0 ), config_fname_( "" ) {
    // init configuration
    config_ = new Config();
    std::cout << "Parsing command line.." << std::endl;
    config_->parseCmdLine( argc, argv );
    if( !config_->needsVersion() && !config_->needsHelp() ) {
        std::cout << "Parsing configuration file.." << std::endl;
        config_->parseFile();
        config_fname_ = config_->optionAsString( "config" );
                
        // bug: somehow the default values of the configuration are not 
        // available before a file has been parsed?
        
        // init globals
        // note: no calls to the random number generator may have been made at
        // this point
        generator.seed( config_->optionAsInt( "init_seed" ) );
        uniform_gen_type aux( generator );
        uniform = aux;
        
        stream_ = new StreamManager( this );
        model_ = new Model( this );
    }
}

fluke::Fluke::~Fluke() {
    if( stream_ != 0 ) delete stream_;
    if( model_ != 0 ) delete model_;
    delete config_;
}

int
fluke::Fluke::run() {
    // run the program swiftly and neatly ;)
    if( config_->needsHelp() ) {
        config_->help( std::cout );
    } else if( config_->needsVersion() ) {
        config_->version( std::cout );
    } else if( config_->needsOverview() ) {
        config_->overview( std::cout );
    } else {
        simulate();
    }
    return 0;
}

void
fluke::Fluke::simulate() {
    // how many times?
    int runs = config_->optionAsInt( "runs" );
    
    // init model
    std::cout << "Building.." << std::endl;
    model_->build();
    
    // and start the simulation random seed generator
    generator.seed( config_->optionAsInt( "random_seed" ) );
    uniform_gen_type aux( generator );
    uniform = aux;

    // one simulation is always run 
    doRun();
    for( int rr = 1; rr != runs; ++rr ) {
        reconfigure( rr );
        doRun();
    }
}

void
fluke::Fluke::doRun() {
    std::cout << "Initializing.." << std::endl;
    model_->initialise();
    
    // write simulation parameters to file just be4 we start running
    boost::filesystem::ofstream *aux = 
        stream_->openOutFileStream( std::string( "simulation.cfg" ),
            std::fstream::out );
    *aux << *config_;
    stream_->closeOutFileStream( aux );
    
    // run simulation
    std::cout << "Running.." << std::endl;
    while( !model_->hasEnded() ) {
        model_->step();
    }
    model_->finish();
}

void
fluke::Fluke::reconfigure( int r ) {
    // sort-of rerunning a modified version of the constructor
    // round up
    config_->reset();

    // assuming the file is present...
    std::string bux( config_fname_ );
    std::string::size_type aux = bux.find_last_of( "." );
    bux.replace( aux - 1, 1, boost::lexical_cast< std::string >( r ) );
    config_->parseFile( bux );
    
#ifdef DEBUG
    cout << "! rebuild" << endl;
#endif
    // new output path
    stream_->createSimulationPath();
    // start with correct population
    model_->rebuild();
    // and start the simulation random seed generator
    generator.seed( config_->optionAsInt( "random_seed" ) );
    uniform_gen_type cux( generator );
    uniform = cux;
}

