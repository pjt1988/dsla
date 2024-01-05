#ifndef TIMING_H_
#define TIMING_H_

#include "enums.h"

#include <chrono>
#include <map>
#include <memory>



namespace DSLA{


    class Timer{
        public:
        Timer();
            ~Timer();
            Timer(const Timer&) = delete;
            Timer(Timer&&) = delete;

            void start();
            void stop(const Timings t, bool reset=false);
            void reset(const Timings t);
            void resetAll();
            void printTimes();

            static std::shared_ptr<Timer> Instance();           


        private:
            std::map<Timings, double> _myTimings;
            std::chrono::time_point<std::chrono::high_resolution_clock> _start;
            std::chrono::time_point<std::chrono::high_resolution_clock> _stop;
            bool _running;
            static std::shared_ptr<Timer> _instance;// = nullptr;         
      
    };

    


}

#endif