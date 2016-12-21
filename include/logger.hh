//
// Some observers of the population.
//
// by Anton Crombach, A.B.M.Crombach@bio.uu.nl
//

#ifndef _FLUKE_LOGGER_H_
#define _FLUKE_LOGGER_H_

#include "defs.hh"
#include "observer.hh"
#include "population.hh"
#include "stream_manager.hh"

namespace fluke {

    class LogCsvMutations : public AsyncLogObserver {
        public:
        LogCsvMutations( std::string, StreamManager * );
        virtual ~LogCsvMutations() {}

        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void writeMutations();
        int classifyDsb( int, int );
        int classifyGeneMutation( int, int );
        
        private:
        long time_;
        uint unknown_;
        std::vector< uint > dsbs_, dsbs_raw_, cps_, rms_;
    };
    
    class LogCsvGenes : public LogObserver {
        public:
        LogCsvGenes( std::string, StreamManager *, long );
        virtual ~LogCsvGenes() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogCsvRates : public LogObserver {
        public:
        LogCsvRates( std::string, StreamManager *, long );
        virtual ~LogCsvRates() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };
    
    class LogCsvPrunedDist : public LogObserver {
        public:
        LogCsvPrunedDist( std::string, StreamManager *, long );
        virtual ~LogCsvPrunedDist() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogCsvDistances : public LogObserver {
        public:
        LogCsvDistances( std::string, StreamManager *, long );
        virtual ~LogCsvDistances() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogCsvScores : public LogObserver {
        public:
        LogCsvScores( std::string, StreamManager *, long );
        virtual ~LogCsvScores() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogCsvEnvironment : public AsyncLogObserver {
        public:
        LogCsvEnvironment( std::string, StreamManager * );
        virtual ~LogCsvEnvironment() {}
        
        virtual void update( Subject * );
        virtual void finalize() {}
    };        

    class LogCsvAncestors : public AsyncLogObserver {
        public:
        LogCsvAncestors( std::string, StreamManager * );
        virtual ~LogCsvAncestors() {}
        
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogXmlAgentTrace : public AsyncLogObserver {
        // The trace file has a certain format. Unfortunately xml is a bit of
        // a deception, so the trace file is in csv format:
        //
        // each line is an agent id-tag: birth <tab> x <tab> y <nl>
        public:
        LogXmlAgentTrace( std::string, std::string, StreamManager * );
        virtual ~LogXmlAgentTrace() {}
        
        virtual void update( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        void writeFooter();
        void nextTarget();
        
        private:
        bool first_;
        std::string dname_;
        AgentTag target_;
        boost::filesystem::ifstream *trace_;
    };
        
    class LogXmlGenomes : public LogObserver {
        public:
        LogXmlGenomes( std::string, StreamManager *, long );
        virtual ~LogXmlGenomes() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    class LogXmlEnvGenomes : public AsyncLogObserver {
        public:
        LogXmlEnvGenomes( std::string, StreamManager * );
        virtual ~LogXmlEnvGenomes() {}
        
        virtual void update( Subject * );
        virtual void finalize()
        { writeFooter(); }
        
        private:
        void writeHeader();
        void writeFooter();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

    class LogPopulationSize : public LogObserver {
        public:
        LogPopulationSize( std::string, StreamManager *, long );
        virtual ~LogPopulationSize() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };
    
    class LogCsvPopulationDistances : public LogObserver {
        public:
        LogCsvPopulationDistances( std::string, StreamManager *, long );
        virtual ~LogCsvPopulationDistances() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
    };

    class LogCsvGrid : public LogObserver {
        public:
        LogCsvGrid( std::string, StreamManager *, long );
        virtual ~LogCsvGrid() {}
        
        virtual void doUpdate( Subject * );
        virtual void finalize() {}
        
        private:
        void writeHeader();
        std::string unique_name( long ) const;
        
        private:
        std::string dname_;
    };

}
#endif

