#include "timing.h"
#include "enums.h"

#include <iostream>

namespace DSLA{
    const std::string enumAsString(const Timings& type){
        using enum Timings;
        switch(type){
        case BCSR_COPY:
            return "BCSR::COPY";
        case BCSR_ADD:
            return "BCSR::ADD";
        case BCSR_SUB:
            return "BCSR::SUB";
        case BCSR_SCALE:
            return "BCSR::SCALE";
        case MISC:
            return "Misc";
        default:
            return "undefined";
        }
    }

    std::shared_ptr<Timer> Timer::_instance = nullptr;

    Timer::Timer(): _myTimings(), _start(std::chrono::high_resolution_clock::now()), _stop(std::chrono::high_resolution_clock::now()), _running(false)
        {}

    std::shared_ptr<Timer> Timer::Instance(){
        if(_instance == nullptr){
            _instance = std::make_shared<Timer>();
        }
        return _instance;
    }

    Timer::~Timer(){
        std::cout << "timer dtor\n";
        printTimes();
        //hacky, but works..
    }

    void Timer::start(){
        _running = true;
        _start = std::chrono::high_resolution_clock::now();
    }

    void Timer::stop(const Timings t, bool reset){
        if(_running){
            _stop = std::chrono::high_resolution_clock::now();
            const auto duration = std::chrono::duration_cast<std::chrono::microseconds> (_stop - _start).count();
            _running = false;
            if(!reset){
                if(_myTimings.count(t) == 0)
                    _myTimings[t] = duration;
                else
                    _myTimings[t] += duration;
            }
            else
                _myTimings[t] = duration;

        }else{
            std::cout << "not running?\n";
            //throw exception or something?
        }
    }

    void Timer::reset(const Timings t){
        _myTimings[t] = 0.0;
    }

    void Timer::resetAll(){
        _myTimings.clear();
    }

    void Timer::printTimes(){
        if(_myTimings.size()){
            for(const auto& it : _myTimings){
                //TODO: replace with std::format
                std::cout << enumAsString(it.first) << " - " << it.second << " us\n";
            }
        }else
            std::cout << "timing is empty...\n";

    }



    



}
