//
// XML sax reader for populations
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "defs.hh"
#include "factory.hh"
#include "population_reader.hh"

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

fluke::PopulationReader::PopulationReader( Factory *ft ) 
    : factory_( ft ), sawErrors_( false ), done_( false ), class_( false ) {
    population_ = std::vector< Agent* >();
    locations_ = std::vector< Location >();
    chr_ = 0;
    chromo_ = 0;
    genome_ = 0;
    agent_ = 0;
    type_ = 0;
    cp_gene_ = -1.0;
    rm_gene_ = -1.0;
    cp_tp_ = -1.0;
    rm_tp_ = -1.0;
    rm_ltr_ = -1.0;
    dsb_ = -1.0;
}

void
fluke::PopulationReader::startElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname,
        const Attributes &attrs ) {
    // now check type of element
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // create tag
        char *bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        // not reading type coz it's overriden in program
        type_ = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        me_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        me_.x = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        me_.y = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 4 ) );
        me_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        chr_ = new std::list< ChromosomeElement* >();
    } else if( aux == "class" ) {
        // we know what kind of agent to create
        class_ = true;
    } else if( aux == "parent" ) {
        // create parent tag
        char* bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        mother_.time = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        mother_.x = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        mother_.y = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );        
        bux = XMLString::transcode( attrs.getValue( 3 ) );
        mother_.i = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );        
    } else if( aux == "dstream" ) {
        // create another downstream
        char* bux = XMLString::transcode(
            attrs.getValue( static_cast< uint >( 0 ) ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        if( attrs.getLength() == 1 ) {
            chr_->push_back( new OrdinaryDownstream( cux ) );
        } else {
            bux = XMLString::transcode( attrs.getValue( 1 ) );
            int eux = boost::lexical_cast< int >( bux );
            XMLString::release( &bux );
            chr_->push_back( new ModuleDownstream( cux, eux ) );
        }
    } else if( aux == "repeat" ) {
        // create repeat
        chr_->push_back( new Repeat() );
    } else if( aux == "tposon" ) {
        // create tposon
        char* bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        int cux = boost::lexical_cast< int >( bux );
        XMLString::release( &bux );
        chr_->push_back( new Retroposon( cux ) );
    } else if( aux == "copy" ) {
        // set copy rates
        char* bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        cp_gene_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        cp_tp_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );                
    } else if( aux == "remove" ) {
        // set removal rates
        char* bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        rm_gene_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 1 ) );
        rm_tp_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );
        bux = XMLString::transcode( attrs.getValue( 2 ) );
        rm_ltr_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );        
    } else if( aux == "break" ) {
        // set break rates
        char* bux = XMLString::transcode( 
            attrs.getValue( static_cast< uint >( 0 ) ) );
        dsb_ = boost::lexical_cast< double >( bux );
        XMLString::release( &bux );
    }
}

void
fluke::PopulationReader::characters( const XMLCh* const chars, 
        const unsigned int length) {
    char *xxx = XMLString::transcode( chars );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( class_ ) {
        agent_class_ = aux;
        class_ = false;
    } 
}

void
fluke::PopulationReader::endElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname ) {
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // ready to create the right agent
        if( agent_class_ == "ModuleAgent" ) {
            //std::cout << "Type: " << type_ << std::endl;
            agent_ = new ModuleAgent( type_, genome_ );
            agent_->myTag( me_ );
            agent_->parentTag( mother_ );
            // and insert in population
            Location bux;
            bux.x = me_.x;
            bux.y = me_.y;
            population_.push_back( agent_ );
            locations_.push_back( bux );
            resetDocument();
        } else {
            throw "Don't know this agent! Peace out..";
        }
    } else if( aux == "genome" ) {
        // create genome (with only one chromosome)
        std::list< Chromosome* > *ll = new std::list< Chromosome* >();
        ll->push_back( chromo_ );
        genome_ = new Genome( ll );
    } else if( aux == "chromosome" ) {
        // create chromosome
        chromo_ = new Chromosome( 0, chr_ );
        chromo_->newBsiteRate( 
            factory_->conf_->optionAsDouble( "new_bsite", type_ ) );
        chromo_->copyBsiteRate( 
            factory_->conf_->optionAsDouble( "cp_bsite", type_ ) );
        chromo_->removeBsiteRate( 
            factory_->conf_->optionAsDouble( "rm_bsite", type_ ) );
        if( cp_gene_ > 0.0 ) {
            chromo_->copyGeneRate( cp_gene_ );
        } else {
            chromo_->copyGeneRate(
                factory_->conf_->optionAsDouble( "cp_gene", type_ ) );
        }
        if( rm_gene_ > 0.0 ) {
            chromo_->removeGeneRate( rm_gene_ );
        } else {
            chromo_->removeGeneRate(
                factory_->conf_->optionAsDouble( "rm_gene", type_ ) );
        }
        if( cp_tp_ > 0.0 ) {
            chromo_->copyRetroposonRate( cp_tp_ );
        } else {
            chromo_->copyRetroposonRate(
                factory_->conf_->optionAsDouble( "cp_tp", type_ ) );
        }
        if( rm_tp_ > 0.0 ) {
            chromo_->removeRetroposonRate( rm_tp_ );
        } else {
            chromo_->removeRetroposonRate(
                factory_->conf_->optionAsDouble( "rm_tp", type_ ) );
        }
        if( dsb_ > 0.0 ) {
            chromo_->recombinationRate( dsb_ );
        } else {
            chromo_->recombinationRate(
                factory_->conf_->optionAsDouble( "dsb_recombination", type_ ) );
        }
        if( rm_ltr_ > 0.0 ) {
            chromo_->removeRepeatRate( rm_ltr_ );
        } else {
            chromo_->removeRepeatRate(
                factory_->conf_->optionAsDouble( "rm_ltr", type_ ) );
        }
        chromo_->newRetroposonRate(
            factory_->conf_->optionAsDouble( "new_tp", type_ ) );
        chromo_->mutationStep( 
            factory_->conf_->optionAsDouble( "mut_step", type_ ) );
        chromo_->dsbStep( 
            factory_->conf_->optionAsDouble( "dsb_step", type_ ) );
        chromo_->retroStep( 
            factory_->conf_->optionAsDouble( "retro_step", type_ ) );
        chromo_->mutationRate( 
            factory_->conf_->optionAsDouble( "mut_rate", type_ ) );
        chromo_->mutationScheme( factory_->mutateRates( type_ ) );
    } else if( aux == "simulation" ) {
        done_ = true;
    }
}

void
fluke::PopulationReader::resetDocument() {
    // just set to zero, we do not own this population
    agent_ = 0;
    genome_ = 0;
    chromo_ = 0;
    chr_ = 0;
    class_ = false;
}

void
fluke::PopulationReader::error( const SAXParseException &e ) {
    sawErrors_ = true;
}

void 
fluke::PopulationReader::fatalError( const SAXParseException &e ) {
    sawErrors_ = true;
}

void
fluke::PopulationReader::warning( const SAXParseException &e ) {}

void
fluke::PopulationReader::resetErrors() {
    sawErrors_ = false;
}

boost::tuple< std::vector< fluke::Agent* >, std::vector< fluke::Location > >
fluke::PopulationReader::population() {
    return boost::make_tuple( population_, locations_ );
}
