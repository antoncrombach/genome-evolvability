//
// XML sax reader for agents
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "defs.hh"
#include "factory.hh"
#include "agent_reader.hh"

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

fluke::AgentReader::AgentReader( Factory *ft ) 
    : factory_( ft ), sawErrors_( false ), done_( false ), class_( false ), 
      ess_( false ), mods_( false ) {
    chr_ = new std::list< ChromosomeElement* >();
    chromo_ = 0;
    genome_ = 0;
    agent_ = 0;
    type_ = 0;
}

void
fluke::AgentReader::startElement( const XMLCh* const uri, 
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
        //type_ = boost::lexical_cast< int >( bux );
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
    }
}

void
fluke::AgentReader::characters( const XMLCh* const chars, 
        const uint length) {
    char *xxx = XMLString::transcode( chars );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( class_ ) {
        agent_class_ = aux;
        class_ = false;
    } 
}

void
fluke::AgentReader::endElement( const XMLCh* const uri, 
        const XMLCh* const localname, const XMLCh* const qname ) {
    char *xxx = XMLString::transcode( localname );
    std::string aux( xxx );
    XMLString::release( &xxx );
    if( aux == "agent" ) {
        // ready to create the right agent
        if( agent_class_ == "ModuleAgent" ) {
            agent_ = new ModuleAgent( type_, genome_ );
            agent_->myTag( me_ );
            agent_->parentTag( mother_ );
        } else {
            throw "Don't know this agent! Peace out..";
        }
        done_ = true;
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
        chromo_->copyGeneRate( 
            factory_->conf_->optionAsDouble( "cp_gene", type_ ) );
        chromo_->removeGeneRate( 
            factory_->conf_->optionAsDouble( "rm_gene", type_ ) );
        chromo_->copyRetroposonRate( 
            factory_->conf_->optionAsDouble( "cp_tp", type_ ) );
        chromo_->removeRetroposonRate( 
            factory_->conf_->optionAsDouble( "rm_tp", type_ ) );
        chromo_->newRetroposonRate(
            factory_->conf_->optionAsDouble( "new_tp", type_ ) );
        chromo_->recombinationRate( 
            factory_->conf_->optionAsDouble( "dsb_recombination", type_ ) );
        chromo_->removeRepeatRate( 
            factory_->conf_->optionAsDouble( "rm_ltr", type_ ) );
        chromo_->mutationScheme( factory_->mutateRates( type_ ) );
    }
}

void
fluke::AgentReader::resetDocument() {
    // just set to zero, we do not own this agent
    chromo_ = 0;
    genome_ = 0;
    agent_ = 0;
}

void
fluke::AgentReader::error( const SAXParseException &e ) {
    sawErrors_ = true;
}

void 
fluke::AgentReader::fatalError( const SAXParseException &e ) {
    sawErrors_ = true;
}

void
fluke::AgentReader::warning( const SAXParseException &e ) {}

void
fluke::AgentReader::resetErrors() {
    sawErrors_ = false;
}

fluke::Agent*
fluke::AgentReader::agent() {
    return agent_;
}
