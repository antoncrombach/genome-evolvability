//
// Implementation of a few observers that dump their data in files.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "logger.hh"
#include "population.hh"
#include "module_agent.hh"
#include "genome.hh"
#include "chromosome.hh"
#include "retroposon.hh"
#include "downstream.hh"
#include "bsite.hh"
#include "environment.hh"
#include "duo_agent.hh"

//
// Counting double strand breaks and other mutations
//
fluke::LogCsvMutations::LogCsvMutations(
    std::string fname, StreamManager *s ) 
    : AsyncLogObserver( s ), time_( -1 ), unknown_( 0 ),
      dsbs_( 6, 0 ), dsbs_raw_( 6, 0 ), cps_( 4, 0 ), rms_( 4, 0 ) {
    openLog( fname );
    *log_ << "# WARNING Only handles one module\n";
    *log_ << "# col 1: time\n";
    *log_ << "# col 2-6: pos/neg, neg/neg, neu/neg, neu/neu, 0/neu unknown\n";
    *log_ << "# col 7-11: idem, raw counts\n";
    *log_ << "# col 12-15: pos, neu, neg (cp)\n";
    *log_ << "# col 16-19: pos, neu, neg (rm)\n";
}

void
fluke::LogCsvMutations::writeMutations() {
    if( std::accumulate( dsbs_.begin(), dsbs_.end(), 0u ) > 0 ||
        std::accumulate( cps_.begin(), cps_.end(), 0u ) > 0 ||
        std::accumulate( rms_.begin(), rms_.end(), 0u ) > 0 ) {
        // new timestep, write old stuff if nonzero
        *log_ << time_ << "\t";
        std::copy( dsbs_.begin(), dsbs_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        std::copy( dsbs_raw_.begin(), dsbs_raw_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        std::copy( cps_.begin(), cps_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        std::copy( rms_.begin(), rms_.end(), 
            std::ostream_iterator< uint >( *log_, "\t" ) );
        *log_ << std::endl;
        
        // do not want accumulative: reset to zero
        std::fill( dsbs_.begin(), dsbs_.end(), 0 );
        std::fill( dsbs_raw_.begin(), dsbs_raw_.end(), 0 );
        std::fill( cps_.begin(), cps_.end(), 0 );
        std::fill( rms_.begin(), rms_.end(), 0 );
    }
}    

void
fluke::LogCsvMutations::update( Subject *s ) {
    // We are receiving two agents..
    DuoAgent *da = dynamic_cast< DuoAgent* >( s );
    ModuleAgent *m1 = dynamic_cast< ModuleAgent* >( da->first );
    ModuleAgent *m2 = dynamic_cast< ModuleAgent* >( da->second );
    if( m1 && m2 ) {
        // get the time
        if( m1->myTag().time > time_ ) {
            writeMutations();
            time_ = m1->myTag().time;
        }
        
        // dsb, copy gene, remove gene
        int da, db, c1, r1;
        // dsbs per chromosome d = dsb
        boost::tie( da, db ) = m1->nrDsbParent();
        // other mutations c = copy, r = remove
        std::vector< uint > mm( m1->nrMutations() );
        c1 = mm[ Chromosome::CP_G ];
        r1 = mm[ Chromosome::RM_G ];
        // length differences
        int ds1 = m1->sizeParent() - m1->size();
        int ds2 = m2->sizeParent() - m2->size();
        // gene dist differences
        int dd1 = m1->distanceParent() - m1->distance();
        int dd2 = m2->distanceParent() - m2->distance();
        
        if( da == 0 || db == 0 ) {
            // indels
            if( c1 == 0 && r1 == 0 ) {
                // nothing happened
            } else if( c1 == 0 && r1 != 0 ) {
                // only rm
                int aux = classifyGeneMutation( dd1, ds1 );
                if( aux != -1 ) {
                    rms_[ aux ] += r1;
                }
                aux = classifyGeneMutation( dd2, ds2 );
                if( aux != -1 ) {
                    rms_[ aux ] += r1;
                }
            } else if( c1 != 0 && r1 == 0 ) {
                // only cp
                int aux = classifyGeneMutation( dd1, ds1 );
                if( aux != -1 ) {
                    cps_[ aux ] += c1;
                }
                aux = classifyGeneMutation( dd2, ds2 );
                if( aux != -1 ) {
                    cps_[ aux ] += c1;
                }
            } else {
                // both aargh
                if( ds1 > 0 ) {
                    // cp?
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        cps_[ aux ] += c1;
                    }
                } else if( ds1 < 0 ) {
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        rms_[ aux ] += r1;
                    }
                }                    
                if( ds2 > 0 ) {
                    // cp?
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        cps_[ aux ] += c1;
                    }
                } else if( ds2 < 0 ) {
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        rms_[ aux ] += r1;
                    }
                }
            }
            // and the dsbs that might have happened are added to 0/neu
            if( da + db != 0 ) ++dsbs_[ 4 ];
            dsbs_raw_[ 4 ] += ( da + db );
        } else {
            if( c1 == 0 && r1 == 0 ) {
                // dsb
                int aux = classifyDsb( dd1, dd2 );
                if( aux != -1 ) {
                    ++dsbs_[ aux ];
                    dsbs_raw_[ aux ] += ( da + db );
                } else {
                    ++unknown_;
                }
            } else {
                // theoretically possible, not seen in sims yet
                // assuming dsb has largest impact, saying cp/rm are neutral
                int aux = classifyDsb( dd1, dd2 ); 
                ++dsbs_[ aux ];
                dsbs_raw_[ aux ] += ( da + db );
                cps_[ 1 ] += c1;
                rms_[ 1 ] += r1;
            }
        } 
 
        /*
        if( c1 + r1 > 0 && da >= 1 && db >= 1 ) {
            // theoretically possible, not seen in sims yet
            // assuming dsb has largest impact, saying cp/rm are neutral
            ++dsbs_[ classifyDsb( dd1, dd2 ) ];
            cps_[ 1 ] += c1;
            rms_[ 1 ] += r1;
        } else if( c1 + r1 > 0 && da + db < 2 ) {
            if( c1 > 0 && r1 == 0 ) {
                int aux = classifyGeneMutation( dd1, ds1 );
                if( aux != -1 ) {
                    cps_[ aux ] += c1;
                }
                aux = classifyGeneMutation( dd2, ds2 );
                if( aux != -1 ) {
                    cps_[ aux ] += c1;
                }
            } else if( c1 == 0 && r1 > 0 ) {
                int aux = classifyGeneMutation( dd1, ds1 );
                if( aux != -1 ) {
                    rms_[ aux ] += r1;
                }
                aux = classifyGeneMutation( dd2, ds2 );
                if( aux != -1 ) {
                    rms_[ aux ] += r1;
                }
            } else { // c1 > 0 && r1 > 0
                // not perfect, but will do...
                if( ds1 > 0 ) {
                    // cp?
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        cps_[ aux ] += c1;
                    }
                } else if( ds1 < 0 ) {
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        rms_[ aux ] += r1;
                    }
                }                    
                if( ds2 > 0 ) {
                    // cp?
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        cps_[ aux ] += c1;
                    }
                } else if( ds2 < 0 ) {
                    int aux = classifyGeneMutation( dd1, ds1 );
                    if( aux != -1 ) {
                        rms_[ aux ] += r1;
                    }
                }                    
            }
        } else if( c1 + r1 == 0 && da >= 1 && db >= 1 ) {
            // only dsb
            int aux = classifyDsb( dd1, dd2 );
            if( aux != -1 ) {
                ++dsbs_[ aux ];
            } else {
                ++unknown_;
            }
        }*/
   }
}

int
fluke::LogCsvMutations::classifyGeneMutation( int dd, int ds ) {
    // only cp xor rm events; pos, neu, neg?
    int result = -1;
    if( ds != 0 ) {
        if( dd > 0 ) {
            result = 0;
        } else if( dd == 0 ) {
            result = 1;
        } else if( dd < 0 ) {
            result = 2;
        }
    }
    return result;
}

int
fluke::LogCsvMutations::classifyDsb( int dd1, int dd2 ) {
    // 4 possibilities, I think
    int result = -1;
    if( dd1 > 0 || dd2 > 0 ) {
        // one pos, other assumed neg
        result = 0;
    } else if( dd1 < 0 && dd2 < 0 ) {
        result = 1;
    } else if( ( dd1 == 0 && dd2 < 0 ) || ( dd1 < 0 && dd2 == 0 ) ) {
        result = 2;
    } else if( dd1 == 0 && dd2 == 0 ) {
        result = 3;
    }
    return result;
}

//
// Simple csv observer
//
fluke::LogCsvGenes::LogCsvGenes( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvGenes::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // first find out dimensions
    uint mod = 0;
    ModuleAgent *ma = dynamic_cast< ModuleAgent * >( am.begin()->first );
    if( ma ) {
        // add 2 for essential genes (index 0) and retroposons (index last)
        mod = ma->nrModules() + 3;
    }
    // init vector
    std::vector< std::vector< uint > > genes( mod );
    
    // fill vector
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        // counting genes in module agents
        ModuleAgent *ma = dynamic_cast< ModuleAgent* >( i->first );
        if( ma ) {
            // first essential genes
            std::vector< uint > aux = ma->nrEssentialGenes();
            genes[ 0 ].push_back( std::accumulate( 
                    aux.begin(), aux.end(), 0 ) );
            // and module genes
            std::vector< std::vector< uint > > bux = ma->nrModuleGenes();
            for( uint ii = 0; ii < bux.size(); ++ii ) {
                genes[ ii + 1 ].push_back( std::accumulate( 
                    bux[ ii ].begin(), bux[ ii ].end(), 0 ) );
            }
            // and add transposons
            int gux = ma->nrRetroposons();
            genes[ bux.size() + 1 ].push_back( ( gux < 0 )? 0: gux );
            // and repeats
            int fux = ma->nrRepeats();
            genes[ bux.size() + 2 ].push_back( ( fux < 0 )? 0: fux );
        }
    }
        
    // calculate max, mean, variance for essential genes, modules and tposons
    typedef  std::vector< std::vector< uint > >::iterator vv_iter;
    for( vv_iter k = genes.begin(); k != genes.end(); ++k ) {
        // calculate some stuff
        double aux = mean( k->begin(), k->end() ); 
        std::stable_sort( k->begin(), k->end() );
        // median
        uint n = k->size();
        double median = 0.0;
        if( n % 2 == 0 ) {
            median = ( ( *k )[ n / 2 ] + ( *k )[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = ( *k )[ ( n / 2 ) ];
        }

        // write max
        *log_ << k->back() << "\t";
        // write mean
        *log_ << aux << "\t";
        // write median
        *log_ << median << "\t";
        *log_ << variance( k->begin(), k->end(), aux ) << "\t";
    }
    *log_ << "\n";
    log_->flush();
}

void
fluke::LogCsvGenes::writeHeader() {
    *log_ << "# max, mean, median, variance per essential and module genes\n";
}

//
// Simple csv observer for tracking the evolving rates
//
fluke::LogCsvRates::LogCsvRates( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvRates::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // fill vector
    std::vector< std::vector< double > > rts( 6, std::vector< double >() );
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        // getting mutation rate values in module agents
        ModuleAgent *ma = dynamic_cast< ModuleAgent* >( i->first );
        if( ma ) {
            // assuming one chromosome...
            const Chromosome *aux = ma->genome().chromosomes().front();
            rts[ 0 ].push_back( aux->copyGeneRate() );
            rts[ 1 ].push_back( aux->removeGeneRate() );
            rts[ 2 ].push_back( aux->recombinationRate() );
            rts[ 3 ].push_back( aux->copyRetroposonRate() );
            rts[ 4 ].push_back( aux->removeRetroposonRate() );
            rts[ 5 ].push_back( aux->removeRepeatRate() );
        }
    }
        
    // calculate mean, variance
    typedef  std::vector< std::vector< double > >::iterator vv_iter;
    for( vv_iter k = rts.begin(); k != rts.end(); ++k ) {
        double aux = mean( k->begin(), k->end() );
        std::stable_sort( k->begin(), k->end() );
        // median
        uint n = k->size();
        double median = 0.0;
        if( n % 2 == 0 ) {
            median = ( ( *k )[ n / 2 ] + ( *k )[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = ( *k )[ ( n / 2 ) ];
        }
        double cux = k->front();
        double dux = k->back();        
        double bux = variance( k->begin(), k->end(), aux );
        *log_ << cux << "\t" << median << "\t" << aux 
              << "\t" << dux << "\t" << bux << "\t";
    }
    *log_ << "\n";
    log_->flush();
}

void
fluke::LogCsvRates::writeHeader() {
    *log_ << "# min, median, avg, max, var for cp gene, rm, dsb\n";
}


//
// Simple csv observer for the pruned distances
//
fluke::LogCsvPrunedDist::LogCsvPrunedDist( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvPrunedDist::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the distances
    int max_dist = ModuleAgent::maxDistance(); 
    std::vector< double > distances;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        double aux = i->first->distance();
        if( aux < max_dist ) {
            distances.push_back( aux );
        }
    }
    
    // calculate min, median, mean and variance
    if( !distances.empty() ) {
        std::stable_sort( distances.begin(), distances.end() ); 
        double min = distances.front();
        double avg = mean( distances.begin(), distances.end() ); 
        
        // median
        double median = 0.0;
        uint n = distances.size();
        if( n % 2 == 0 ) {
            median = ( distances[ n / 2 ] + distances[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = distances[ ( n / 2 ) ];
        }
        
        // std dev
        double aux = 0.0;
        typedef std::vector< double >::iterator vit;
        for( vit i = distances.begin(); i != distances.end(); ++i ) {
            aux += ( *i - avg ) * ( *i - avg );
        }
        double sdev = sqrt( aux / ( n - 1 ) );
        
        *log_ << min << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
        log_->flush();
    } else {
        *log_ << "# Empty grid...\n";
        log_->flush();
    }
}

void
fluke::LogCsvPrunedDist::writeHeader() {
    *log_ << "# min, median, mean, variance of genotypic distances\n";
}

//
// Simple csv observer for the distances
//
fluke::LogCsvDistances::LogCsvDistances( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the distances
    std::vector< double > distances;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        distances.push_back( i->first->distance() );
    }
    
    // calculate min, median, mean and variance
    if( !distances.empty() ) {
        std::stable_sort( distances.begin(), distances.end() ); 
        double min = distances.front();
        double avg = mean( distances.begin(), distances.end() ); 
        
        // median
        double median = 0.0;
        uint n = distances.size();
        if( n % 2 == 0 ) {
            median = ( distances[ n / 2 ] + distances[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = distances[ ( n / 2 ) ];
        }
        
        // std dev
        double aux = 0.0;
        typedef std::vector< double >::iterator vit;
        for( vit i = distances.begin(); i != distances.end(); ++i ) {
            aux += ( *i - avg ) * ( *i - avg );
        }
        double sdev = sqrt( aux / ( n - 1 ) );
        
        *log_ << min << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
        log_->flush();
    } else {
        *log_ << "# Empty grid...\n";
        log_->flush();
    }
}

void
fluke::LogCsvDistances::writeHeader() {
    *log_ << "# min, median, mean, variance of genotypic distances\n";
}


//
// Simple csv observer for the scores
//
fluke::LogCsvScores::LogCsvScores( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvScores::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::vector< double > scores;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        scores.push_back( i->first->score() );
    }
    
    // calculate max, median, mean and variance
    if( !scores.empty() ) {
        std::stable_sort( scores.begin(), scores.end() ); 
        double max = scores.back();
        double avg = mean( scores.begin(), scores.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores.size();
        if( n % 2 == 0 ) {
            median = ( scores[ n / 2 ] + scores[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores[ ( n / 2 ) ];
        }
        
        // std dev
        double aux = 0.0;
        typedef std::vector< double >::iterator vit;
        for( vit i = scores.begin(); i != scores.end(); ++i ) {
            aux += ( *i - avg ) * ( *i - avg );
        }
        double sdev = sqrt( aux / ( n - 1 ) );
        
        *log_ << max << "\t" << median << "\t" << avg << "\t" << sdev << "\n";
        log_->flush();
    } else {
        *log_ << "# Empty grid...\n";
        log_->flush();
    }
}

void
fluke::LogCsvScores::writeHeader() {
    *log_ << "# max, median, mean, variance of raw fitness scores\n";
}


//
// Simple csv observer for the environment
//
fluke::LogCsvEnvironment::LogCsvEnvironment( 
        std::string fname, StreamManager *s ) : AsyncLogObserver( s ) {
    openLog( fname );
    *log_ << "# time, new environment\n";
}

void
fluke::LogCsvEnvironment::update( Subject *s ) {
    // Can this be done more elegantly?
    Environment *env = dynamic_cast< ConstantEnvironment * >( s );
    if( !env ) env = dynamic_cast< PoissonEnvironment * >( s );
    if( !env ) env = dynamic_cast< PeriodicEnvironment * >( s );
    // log change from which to which state
    *log_ << env->model()->now() << "\t" << env->expectedCopies( 0 ) << "\t"
          << env->expectedCopies( 1 ) << "\n";
    log_->flush();
}


//
// Simple csv observer for ancestor tracing
//
fluke::LogCsvAncestors::LogCsvAncestors( 
    std::string fname, StreamManager *s ) : AsyncLogObserver( s ) {
    openLog( fname );
    writeHeader();
}

void 
fluke::LogCsvAncestors::update( Subject *s ) {
    Agent *ag = dynamic_cast< Agent * >( s );
    AgentTag child = ag->myTag();
    AgentTag mother = ag->parentTag();
    // log birth of new agent
    *log_ << child.time << "\t" << child.x << "\t" << child.y
        << "\t" << child.i << "\t-> ";
    *log_ << mother.time  << "\t" << mother.x << "\t" << mother.y 
        << "\t" << mother.i << "\n";
}

void 
fluke::LogCsvAncestors::writeHeader() {
    *log_ << "# child ( birth, x, y ) -> parent ( birth, x, y )\n";
}


//
// After one run, we can trace back agents to the beginning and then follow
// their development
//
fluke::LogXmlAgentTrace::LogXmlAgentTrace( std::string dname, 
    std::string tracefile, StreamManager *s ) 
    : AsyncLogObserver( s ) {
    first_ = true;
    dname_ = dname;
    s->openPath( dname_ );
    trace_ = s->openInFileStream( tracefile, std::fstream::in );
    nextTarget();
}

void
fluke::LogXmlAgentTrace::update( Subject *s ) {
    ModuleAgent *ag = dynamic_cast< ModuleAgent * >( s );
    AgentTag tag = ag->myTag();
    if( tag == target_ ) {
        // devise name
        std::string fname = tag.str();
        // write genome
        openLog( dname_ + "/" + fname + ".xml" );
        writeHeader();
        int aa, ab;
        boost::tie( aa, ab ) = ag->nrDsbParent();
        std::vector< uint > bb( ag->nrMutations() );
        *log_ << "<mutations dsb_a=\"" << aa
            << "\" dsb_b=\"" << ab 
            << "\" cp_g=\"" << bb[ Chromosome::CP_G ] 
            << "\" rm_g=\"" << bb[ Chromosome::RM_G ] << "\"/>\n";
        *log_ << *ag;
        writeFooter();
        closeLog();
        // and move on
        nextTarget();
    }
}

void
fluke::LogXmlAgentTrace::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
fluke::LogXmlAgentTrace::writeFooter() {
    *log_ << "</simulation>\n";
}

void
fluke::LogXmlAgentTrace::nextTarget() {
    // checks for input / output necessary?
    *trace_ >> target_.time >> target_.x >> target_.y >> target_.i;
}


//
// Another class
// 
fluke::LogXmlGenomes::LogXmlGenomes(
        std::string dname, StreamManager *s, long i ) : LogObserver( s, i ) {
    // create dir with given name
    // within dir, use some logical name
    // f.i. timestep_xcoordinate_ycoordinate.dot
    dname_ = dname;
    s->openPath( dname_ );
}

void
fluke::LogXmlGenomes::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n";

    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        // output the genome agents
        *log_ << *( i->first );
    }
    // end of file
    *log_ << "</generation>\n";
    closeLog();
}

void
fluke::LogXmlGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
fluke::LogXmlGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
fluke::LogXmlGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}

//
// Another class, yet asynchronous
// 
fluke::LogXmlEnvGenomes::LogXmlEnvGenomes( std::string dname, StreamManager *s ) 
    : AsyncLogObserver( s ) {
    // create dir with given name
    // within dir, use some logical name
    // f.i. timestep_xcoordinate_ycoordinate.dot
    dname_ = dname;
    s->openPath( dname_ );
}

void
fluke::LogXmlEnvGenomes::update( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    *log_ << "<generation time=\"" << pop->generation() << "\">\n";

    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        // output the genome agents
        *log_ << *( i->first );
    }
    // end of file
    *log_ << "</generation>\n";
    //writeFooter();
    closeLog();
}

void
fluke::LogXmlEnvGenomes::writeHeader() {
    *log_ << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
          << "<simulation fluke_version=\"" << VERSION << "\">\n";
}

void
fluke::LogXmlEnvGenomes::writeFooter() {
    *log_ << "</simulation>\n";
}

std::string
fluke::LogXmlEnvGenomes::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".xml";
    return result.str();
}


//
// Yet another class
//
fluke::LogPopulationSize::LogPopulationSize(
    std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogPopulationSize::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::map< int, int > subpop;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        int aux = i->first->type();
        ++subpop[ aux /* + 1 */ ];
    }
    
    // and write to log
    /*for( std::map< int, int >::iterator i = subpop.begin();
        i != subpop.end(); ++i ) {
        *log_ << i-> first << " " << i->second << "\t";
    }*/
    *log_ << "1 " << subpop[ 1 ] << "\t";
    *log_ << "2 " << subpop[ 2 ];
    *log_ << std::endl;
}

void
fluke::LogPopulationSize::writeHeader() {
    *log_ << "# population size per agent type\n";
}

//
// Simple csv observer for the population scores
//
fluke::LogCsvPopulationDistances::LogCsvPopulationDistances( 
        std::string fname, StreamManager *s, long i ) : LogObserver( s, i ) {
    openLog( fname );
    writeHeader();
}

void
fluke::LogCsvPopulationDistances::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_map am = pop->map();

    // get all the scores
    std::vector< double > scores_one;
    std::vector< double > scores_two;
    for( Population::map_ag_iter i = am.begin(); i != am.end(); ++i ) {
        if( i->first->type() == 1 ) {
            scores_one.push_back( i->first->distance() );
        } else if( i->first->type() == 2 ) {
            scores_two.push_back( i->first->distance() );
        }
    }
    
    // calculate min, median, mean and variance for score_one
    if( !scores_one.empty() ) {
        std::stable_sort( scores_one.begin(), scores_one.end() ); 
        double min = scores_one.front();
        double avg = mean( scores_one.begin(), scores_one.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores_one.size();
        if( n % 2 == 0 ) {
            median = ( scores_one[ n / 2 ] + scores_one[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores_one[ ( n / 2 ) ];
        }
        
        // std dev
        double aux = 0.0;
        typedef std::vector< double >::iterator vit;
        for( vit i = scores_one.begin(); i != scores_one.end(); ++i ) {
            aux += ( *i - avg ) * ( *i - avg );
        }
        double sdev = sqrt( aux / ( n - 1 ) );
        
        *log_ << "1\t" << min << "\t" << median 
              << "\t" << avg << "\t" << sdev << "\t";
    } else {
        *log_ << "1\t0\t0\t0\t0\t";
    }
    // calculate min, median, mean and variance
    if( !scores_two.empty() ) {
        std::stable_sort( scores_two.begin(), scores_two.end() ); 
        double min = scores_two.front();
        double avg = mean( scores_two.begin(), scores_two.end() ); 
        
        // median
        double median = 0.0;
        uint n = scores_two.size();
        if( n % 2 == 0 ) {
            median = ( scores_two[ n / 2 ] + scores_two[ ( n / 2 ) - 1 ] ) / 2.0;
        } else {
            median = scores_two[ ( n / 2 ) ];
        }
        
        // std dev
        double aux = 0.0;
        typedef std::vector< double >::iterator vit;
        for( vit i = scores_two.begin(); i != scores_two.end(); ++i ) {
            aux += ( *i - avg ) * ( *i - avg );
        }
        double sdev = sqrt( aux / ( n - 1 ) );
        
        *log_ << "2\t" << min << "\t" << median 
              << "\t" << avg << "\t" << sdev << "\n";
        log_->flush();
    } else {
        *log_ << "\n";
    }
}

void
fluke::LogCsvPopulationDistances::writeHeader() {
    *log_ << "# min, median, mean, variance of raw scores per agent type\n";
}

//
// And now a population dumper, at least it dumps a few specific feats
//
fluke::LogCsvGrid::LogCsvGrid( 
        std::string dname, StreamManager *s, long i ) : LogObserver( s, i ) {
    dname_ = dname;
    s->openPath( dname_ );
}

void
fluke::LogCsvGrid::doUpdate( Subject *s ) {
    Population *pop = static_cast< Population * >( s );
    Population::agents_grid grid = pop->grid();
    
    // begin of file
    openLog( unique_name( pop->generation() ) );
    writeHeader();
    // and a label of the feat
    *log_ << "gene distance" << "\t";
    // with the dimensions of the matrix to follow
    int n = grid.shape()[ 0 ];
    int m = grid.shape()[ 1 ];
    *log_ << n << "\t" << m << "\n";
    for( int i = 0; i != n; ++i ) {
        for( int j = 0; j != m; ++j ) {
            ModuleAgent *ag = dynamic_cast< ModuleAgent * >( grid[ i ][ j ] );
            if( ag ) {
                *log_ << std::setprecision( 5 ) << ag->distance() << "\t";
            } else {
                *log_ << std::setprecision( 5 ) << -1 << "\t";
            }
        }
        *log_ << "\n";
    }
    closeLog();
}

void
fluke::LogCsvGrid::writeHeader() {
    *log_ << "# grid size and features in matrix format,";
    *log_ << " separate by newlines\n";
}

std::string
fluke::LogCsvGrid::unique_name( long time ) const {
    std::stringstream result;
    
    result << dname_ << "/";
    // fix: using magic numbers!
    result << "t" << std::setw( 8 ) << std::setfill( '0' ) << time << ".csv";
    return result.str();
}

