//
// Manager for file streams.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "stream_manager.hh"
#include "fluke.hh"
#include "config.hh"


fluke::StreamManager::StreamManager( Fluke *fl ) 
    : fluke_( fl ), basepath_(), simulation_folder_(), 
      in_files_(), out_files_(), open_files_() {}

fluke::StreamManager::~StreamManager() {
    closeAllFileStreams();
    in_files_.clear();
    out_files_.clear();
    open_files_.clear();
}

boost::filesystem::ofstream* 
fluke::StreamManager::openOutFileStream( std::string s, 
        std::ios_base::openmode m = std::fstream::out ) {
    // no error checking yet
    // is this gonna work?
    simulationPath();
    boost::filesystem::ofstream *aux = new boost::filesystem::ofstream();
    /* aux->open( simulation_folder_ / s, m ); */
    
    // patch for >2Gb files
    std::string f_mode;
    if( m & std::fstream::out ) {
        f_mode = "w";
    } else if( m & ( std::fstream::out | std::fstream::app ) ) {
        f_mode = "a";
    } else {
        throw( "Unknown file openmode" );
    }
    FILE *fd = fopen( ( simulation_folder_ / s ).string().c_str(),
        f_mode.c_str() );
    if( fd != NULL ) {
        __gnu_cxx::stdio_filebuf< char > *filebuf = 
            new __gnu_cxx::stdio_filebuf< char >( fd, m ); 
        aux->std::ios::rdbuf( filebuf );
    } else {
        throw "Cannot open file for writing.";
    }
     // end of patch
    
    out_files_.push_back( aux );
    open_files_[ aux ] = fd;
    return aux;
}

boost::filesystem::ifstream*
fluke::StreamManager::openInFileStream( std::string s,
        std::ios_base::openmode m = std::fstream::in ) {
    // no error checking
    // taking just the string as the location, no prepending of a path
    boost::filesystem::ifstream *aux = new boost::filesystem::ifstream();
    aux->open( s, m );
    in_files_.push_back( aux );
    return aux;
}

void 
fluke::StreamManager::closeInFileStream( boost::filesystem::ifstream *fs ) {
    infile_iter aux = std::find( in_files_.begin(), in_files_.end(), fs );
    ( **aux ).close();
    delete *aux;
    in_files_.erase( aux );
}

void 
fluke::StreamManager::closeOutFileStream( boost::filesystem::ofstream *fs ) {
    outfile_iter aux = std::find( out_files_.begin(), out_files_.end(), fs );
    // due to >2Gb patch no automatic flushing of the file?
    ( **aux ).flush();
    delete ( **aux ).std::ios::rdbuf();
    // close the FILE*
    fclose( open_files_[ *aux ] );
    open_files_.erase( *aux );
    ( **aux ).close();
    delete *aux;   
    out_files_.erase( aux );
}

void 
fluke::StreamManager::closeAllFileStreams() {
    for( outfile_iter i = out_files_.begin(); i != out_files_.end(); ++i ) {
        // >2Gb patch
        ( **i ).flush();
        delete ( **i ).std::ios::rdbuf();
        // close the FILE*
        fclose( open_files_[ *i ] );
        open_files_.erase( *i );
        ( **i ).close();
        delete *i;
    }
    out_files_.clear();
    for( infile_iter i = in_files_.begin(); i != in_files_.end(); ++i ) {
        ( **i ).close();
        // and delete everything
        delete *i;
    }    
    in_files_.clear();
}

void
fluke::StreamManager::openPath( std::string path ) {
    // opens path, erases it if existing?
    simulationPath();
    boost::filesystem::create_directory( simulation_folder_ / path );
}

void
fluke::StreamManager::createSimulationPath() {
    basepath_ = "";
    simulation_folder_ = "";
    simulationPath();
}

void 
fluke::StreamManager::simulationPath() {
    using boost::filesystem::directory_iterator;
    using boost::filesystem::create_directory;
    using boost::filesystem::remove;
    using boost::algorithm::split;
    using boost::algorithm::erase_all;

    if( basepath_.empty() and simulation_folder_.empty() ) {
        basepath_ = fluke_->configuration().optionAsString( "log_path" );

        if( symbolic_link_exists( basepath_ / "latest" ) ) {
            remove( basepath_ / "latest" );
        }
        
        int mm = -1;
        directory_iterator end;
        for( directory_iterator i( basepath_ ); i != end; ++i ) {
            if( isDataFolder( i ) ) {
                // only dirs left in format xxx-xxx-xxx-xxx
                std::string bux = i->leaf();
                erase_all( bux, "-" );
                // to int
                int cux = boost::lexical_cast< int >( bux );
                if( cux > mm ) {
                    mm = cux;
                }
            }
        }

        // create new path
        ++mm;
        simulation_folder_ = basepath_ / formatSimFolder( mm );
#ifdef DEBUG
        cout << "! log dir is: " << simulation_folder_.string() << endl;
#endif
        create_directory( simulation_folder_ );

        // and a symbolic link
        symlink( simulation_folder_.leaf().c_str(), 
                 ( basepath_ / "latest" ).string().c_str() );
    }
}

bool
fluke::StreamManager::isDataFolder( 
    boost::filesystem::directory_iterator d ) const {
    // check if folder has format xxx-xxx-xxx-xxx
    using boost::filesystem::directory_iterator;
    using boost::filesystem::symbolic_link_exists;
    using boost::filesystem::is_directory;

    boost::regex expr0( "\\d\\d\\d-\\d\\d\\d-\\d\\d\\d-\\d\\d\\d" );
    boost::regex expr1( "\\d\\d\\d-\\d\\d\\d" );
    
    bool result = false;
    if( is_directory( *d ) && !symbolic_link_exists( *d ) ) {
        std::string bux = d->leaf();
        result = boost::regex_match( bux, expr0 ) || 
            boost::regex_match( bux, expr1 );
    }
    return result;
}

boost::filesystem::path
fluke::StreamManager::formatSimFolder( int n ) {
    // using magic numbers!!
    int k = 6;
    std::stringstream aux;
    aux << std::setfill( '0' ) << std::setw( k ) << n;
    std::string result = aux.str();

    for( int i = k - 3; i > 1; i -= 3 ) {
        result.insert( i, "-" );
    }

    return boost::filesystem::path( result );
}

