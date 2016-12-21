//
// Implementation of a few methods in the subject/observer pattern.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#include "observer.hh"
#include "stream_manager.hh"


fluke::LogObserver::LogObserver( StreamManager *sm, long i )
    : Observer(), stream_manager_( sm ), interval_( i ), val_( 1 ) {
    log_ = 0;
}

fluke::LogObserver::LogObserver( const LogObserver &lo ) 
    : Observer( lo ) {
    stream_manager_ = lo.stream_manager_;
    log_ = lo.log_;
    interval_ = lo.interval_;
    val_ = lo.val_;
}

fluke::LogObserver::~LogObserver() {
    closeLog();
}

void 
fluke::LogObserver::openLog( std::string fname ) {
    log_ = stream_manager_->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
fluke::LogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        stream_manager_->closeOutFileStream( log_ );
        log_ = 0;
    }
}

void
fluke::LogObserver::update( Subject *s ) {
    // template method
    if( val_ == 1 ) {
        val_ = interval_;
        doUpdate( s );
    } else {
        --val_;
    }
}

//
// async log observer
//
fluke::AsyncLogObserver::AsyncLogObserver( StreamManager *sm )
    : Observer(), stream_manager_( sm ) {
    log_ = 0;
}

fluke::AsyncLogObserver::AsyncLogObserver( const AsyncLogObserver &lo ) 
    : Observer( lo ) {
    stream_manager_ = lo.stream_manager_;
    log_ = lo.log_;
}

fluke::AsyncLogObserver::~AsyncLogObserver() {
    closeLog();
}

void 
fluke::AsyncLogObserver::openLog( std::string fname ) {
    log_ = stream_manager_->openOutFileStream( fname, 
            std::fstream::out | std::fstream::app );
}

void 
fluke::AsyncLogObserver::closeLog() {
    if( log_ != 0 ) {
        finalize();
        stream_manager_->closeOutFileStream( log_ );
        log_ = 0;
    }
}

